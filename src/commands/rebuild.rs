use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::fmt::Write as FmtWrite;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::{self, Command};
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use crate::commands::shared::{
    print_contract, resolve_gfa_editor_cli, run_gfa_editor_image, CommandContract, GfaImageExport,
};
use crate::error::OrgraftError;
use crate::commands::resolve::sha256_hex;

mod overlap_only {
    //! Final-stage GFA overlap repair. Node sequences, topology, paths, and FASTA
    //! are immutable here; only an unambiguous L-field overlap can be changed.

    use std::collections::BTreeMap;
    use std::fs;
    use std::path::Path;

    use crate::commands::resolve::sha256_hex;
    use crate::error::OrgraftError;

    use super::{read_fasta, reverse_complement_iupac, Gfa};

    const MAX_NODES: usize = 64;
    const MAX_LINKS: usize = 64;
    const MAX_CANDIDATES_PER_LINK: usize = 4_096;
    const MAX_STATES: usize = 200_000;
    const MAX_TRANSITIONS: usize = 2_000_000;
    const MAX_STARTS: usize = 2_000;
    const MAX_TARGET_BP: usize = 5_000_000;
    const MAX_GFA_BYTES: usize = 100_000_000;

    fn invalid(message: impl Into<String>) -> OrgraftError {
        OrgraftError::InvalidArgument(format!("rebuild graph check: {}", message.into()))
    }

    #[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
    struct Step {
        node: usize,
        reverse: bool,
    }

    impl Step {
        fn flipped(self) -> Self {
            Self { node: self.node, reverse: !self.reverse }
        }

        fn index(self) -> usize {
            self.node * 2 + usize::from(self.reverse)
        }
    }

    #[derive(Clone)]
    struct Edge {
        from: Step,
        to: Step,
        old: usize,
        candidates: Vec<usize>,
        line_index: usize,
    }

    struct PathRecord {
        steps: Vec<Step>,
        overlaps: Vec<Option<usize>>,
    }

    struct Graph {
        names: Vec<String>,
        sequences: Vec<[Vec<u8>; 2]>,
        edges: Vec<Edge>,
        paths: Vec<PathRecord>,
    }

    #[derive(Clone, Copy)]
    struct Arc {
        to: Step,
        edge: usize,
        overlap: usize,
    }

    struct Solution {
        overlaps: Vec<usize>,
        walk: Vec<Step>,
        walk_edges: Vec<usize>,
    }

    struct Proof {
        solution: Solution,
        states: usize,
        transitions: usize,
        starts: usize,
        candidate_counts: Vec<usize>,
    }

    pub(super) struct ReplayProof {
        pub(super) status: &'static str,
        pub(super) path_type: &'static str,
        pub(super) used_segments: usize,
        pub(super) used_links: usize,
        pub(super) states: usize,
        pub(super) transitions: usize,
        pub(super) starts: usize,
    }

    fn parse_overlap(token: &str) -> Result<usize, OrgraftError> {
        let digits = token.strip_suffix('M').ok_or_else(|| invalid("only integer nM overlaps are supported"))?;
        if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
            return Err(invalid("only integer nM overlaps are supported"));
        }
        digits.parse::<usize>().map_err(|_| invalid("overlap length is out of range"))
    }

    fn parse_orientation(token: &str) -> Result<bool, OrgraftError> {
        match token {
            "+" => Ok(false),
            "-" => Ok(true),
            _ => Err(invalid("invalid GFA orientation")),
        }
    }

    fn parse_step(token: &str, names: &BTreeMap<String, usize>) -> Result<Step, OrgraftError> {
        let name = token.get(..token.len().saturating_sub(1))
            .ok_or_else(|| invalid("invalid P step"))?;
        let orientation = token.get(token.len().saturating_sub(1)..)
            .ok_or_else(|| invalid("invalid P step"))?;
        let node = names.get(name).copied().ok_or_else(|| invalid("P references an absent S node"))?;
        Ok(Step { node, reverse: parse_orientation(orientation)? })
    }

    fn body(line: &str) -> &str {
        line.strip_suffix("\r\n").or_else(|| line.strip_suffix('\n')).unwrap_or(line)
    }

    fn exact_candidates(left: &[u8], right: &[u8]) -> Result<Vec<usize>, OrgraftError> {
        // KMP border chain: all exact suffix(left) = prefix(right) lengths, not
        // merely the longest match. NUL cannot occur in validated DNA.
        let mut joined = Vec::with_capacity(right.len() + left.len() + 1);
        joined.extend_from_slice(right);
        joined.push(0);
        joined.extend_from_slice(left);
        let mut pi = vec![0usize; joined.len()];
        for i in 1..joined.len() {
            let mut j = pi[i - 1];
            while j > 0 && joined[i] != joined[j] {
                j = pi[j - 1];
            }
            if joined[i] == joined[j] {
                j += 1;
            }
            pi[i] = j;
        }
        let mut lengths = vec![0];
        let mut n = *pi.last().unwrap_or(&0);
        while n > 0 {
            lengths.push(n);
            if lengths.len() > MAX_CANDIDATES_PER_LINK {
                return Err(invalid("too many exact overlaps for one L edge; human review required"));
            }
            n = pi[n - 1];
        }
        lengths.sort_unstable();
        Ok(lengths)
    }

    fn parse_graph(text: &str, enumerate_overlaps: bool) -> Result<Graph, OrgraftError> {
        let mut raw_sequences = BTreeMap::<String, Vec<u8>>::new();
        let mut raw_edges = Vec::<(String, bool, String, bool, usize, usize)>::new();
        let mut raw_paths = Vec::<(Vec<String>, Vec<Option<usize>>)>::new();
        for (line_index, line) in text.split_inclusive('\n').enumerate() {
            let record = body(line);
            if record.is_empty() {
                continue;
            }
            let fields = record.split('\t').collect::<Vec<_>>();
            match fields[0] {
                "H" => {}
                "S" => {
                    if fields.len() < 3 || fields[1].is_empty() || fields[2].is_empty() || fields[2] == "*" {
                        return Err(invalid("invalid S record"));
                    }
                    let sequence = fields[2].to_ascii_uppercase().into_bytes();
                    if !sequence.iter().all(|b| b"ACGTRYSWKMBDHVN".contains(b)) {
                        return Err(invalid("S node contains non-IUPAC DNA"));
                    }
                    if raw_sequences.insert(fields[1].to_string(), sequence).is_some() {
                        return Err(invalid("duplicate S name"));
                    }
                }
                "L" => {
                    if fields.len() < 6 {
                        return Err(invalid("invalid L record"));
                    }
                    raw_edges.push((
                        fields[1].to_string(), parse_orientation(fields[2])?,
                        fields[3].to_string(), parse_orientation(fields[4])?,
                        parse_overlap(fields[5])?, line_index,
                    ));
                }
                "P" => {
                    if fields.len() < 4 || fields[2].is_empty() || fields[2] == "*" {
                        return Err(invalid("invalid P record"));
                    }
                    let steps = fields[2].split(',').map(str::to_string).collect::<Vec<_>>();
                    let overlaps = if fields[3] == "*" {
                        vec![None; steps.len().saturating_sub(1)]
                    } else {
                        let tokens = fields[3].split(',').collect::<Vec<_>>();
                        if tokens.len() != steps.len().saturating_sub(1) {
                            return Err(invalid("P overlap count differs from its steps"));
                        }
                        tokens.into_iter().map(|token| {
                            if token == "*" { Ok(None) } else { parse_overlap(token).map(Some) }
                        }).collect::<Result<Vec<_>, _>>()?
                    };
                    raw_paths.push((steps, overlaps));
                }
                _ => return Err(invalid("unsupported GFA record may carry an unchecked path")),
            }
        }
        if raw_sequences.is_empty() || raw_edges.is_empty() {
            return Err(invalid("the final circular GFA requires S and L records"));
        }
        if raw_sequences.len() > MAX_NODES || raw_edges.len() > MAX_LINKS {
            return Err(invalid("graph exceeds bounded overlap-only search size"));
        }
        let names = raw_sequences.keys().cloned().collect::<Vec<_>>();
        let indices = names.iter().enumerate().map(|(i, name)| (name.clone(), i))
            .collect::<BTreeMap<_, _>>();
        let sequences = names.iter().map(|name| {
            let forward = raw_sequences[name].clone();
            let reverse = reverse_complement_iupac(
                std::str::from_utf8(&forward).expect("validated ASCII DNA"),
            ).into_bytes();
            [forward, reverse]
        }).collect::<Vec<_>>();
        let mut edges = Vec::new();
        for (from_name, from_reverse, to_name, to_reverse, old, line_index) in raw_edges {
            let from = Step {
                node: indices.get(&from_name).copied().ok_or_else(|| invalid("L references an absent S node"))?,
                reverse: from_reverse,
            };
            let to = Step {
                node: indices.get(&to_name).copied().ok_or_else(|| invalid("L references an absent S node"))?,
                reverse: to_reverse,
            };
            let left = &sequences[from.node][usize::from(from.reverse)];
            let right = &sequences[to.node][usize::from(to.reverse)];
            let candidates = if enumerate_overlaps {
                exact_candidates(left, right)?
            } else {
                if old > left.len() || old > right.len()
                    || (old > 0 && left[left.len() - old..] != right[..old]) {
                    return Err(invalid("an existing L overlap exceeds a segment or has disagreeing bases"));
                }
                vec![old]
            };
            edges.push(Edge { from, to, old, candidates, line_index });
        }
        let paths = raw_paths.into_iter().map(|(steps, overlaps)| {
            Ok(PathRecord {
                steps: steps.iter().map(|token| parse_step(token, &indices))
                    .collect::<Result<Vec<_>, _>>()?,
                overlaps,
            })
        }).collect::<Result<Vec<_>, OrgraftError>>()?;
        Ok(Graph { names, sequences, edges, paths })
    }

    fn canonical_cycle(nodes: &[Step], edges: &[usize]) -> Vec<(usize, bool, usize)> {
        let n = nodes.len();
        let mut best: Option<Vec<(usize, bool, usize)>> = None;
        for offset in 0..n {
            let forward = (0..n).map(|i| {
                let at = (offset + i) % n;
                (nodes[at].node, nodes[at].reverse, edges[at])
            }).collect::<Vec<_>>();
            let reverse = (0..n).map(|i| {
                let at = (offset + n - i) % n;
                (nodes[at].node, !nodes[at].reverse, edges[(at + n - 1) % n])
            }).collect::<Vec<_>>();
            for candidate in [forward, reverse] {
                if best.as_ref().is_none_or(|current| candidate < *current) {
                    best = Some(candidate);
                }
            }
        }
        best.unwrap_or_default()
    }

    struct Search<'a> {
        graph: &'a Graph,
        target: &'a [u8],
        repeated_target: Vec<u8>,
        adjacency: Vec<Vec<Arc>>,
        assignments: Vec<Option<usize>>,
        walk: Vec<Step>,
        walk_edges: Vec<usize>,
        solutions: Vec<Solution>,
        identities: Vec<(Vec<usize>, Vec<(usize, bool, usize)>)>,
        states: usize,
        transitions: usize,
        starts: usize,
        max_overlap: usize,
    }

    impl<'a> Search<'a> {
        fn new(graph: &'a Graph, target: &'a [u8]) -> Self {
            let mut adjacency = vec![Vec::new(); graph.names.len() * 2];
            let mut max_overlap = 0;
            for (edge_index, edge) in graph.edges.iter().enumerate() {
                for &overlap in &edge.candidates {
                    max_overlap = max_overlap.max(overlap);
                    adjacency[edge.from.index()].push(Arc { to: edge.to, edge: edge_index, overlap });
                    adjacency[edge.to.flipped().index()].push(Arc {
                        to: edge.from.flipped(), edge: edge_index, overlap,
                    });
                }
            }
            let mut repeated_target = Vec::with_capacity(target.len() * 3);
            for _ in 0..3 { repeated_target.extend_from_slice(target); }
            Self {
                graph, target, repeated_target, adjacency,
                assignments: vec![None; graph.edges.len()],
                walk: Vec::new(), walk_edges: Vec::new(),
                solutions: Vec::new(), identities: Vec::new(),
                states: 0, transitions: 0, starts: 0, max_overlap,
            }
        }

        fn p_valid(&self, overlaps: &[usize]) -> bool {
            self.graph.paths.iter().all(|path| {
                path.steps.windows(2).enumerate().all(|(i, pair)| {
                    self.graph.edges.iter().enumerate().any(|(edge_index, edge)| {
                        let connected = (edge.from == pair[0] && edge.to == pair[1])
                            || (edge.to.flipped() == pair[0] && edge.from.flipped() == pair[1]);
                        connected && path.overlaps[i].is_none_or(|n| n == overlaps[edge_index])
                    })
                })
            })
        }

        // This is an existence check for an already constructed graph. Alternative
        // repeat links are legal and need not all belong to one accepted circle.
        fn visit_existing(
            &mut self, start: Step, here: Step, position: usize, length: usize,
            used: u64, seen: u64, all_nodes: u64,
        ) -> Result<Option<ReplayProof>, OrgraftError> {
            self.states += 1;
            if self.states > MAX_STATES {
                return Err(invalid("existing-circle search limit reached; human review required"));
            }
            for i in 0..self.adjacency[here.index()].len() {
                self.transitions += 1;
                if self.transitions > MAX_TRANSITIONS {
                    return Err(invalid("existing-circle transition limit reached; human review required"));
                }
                let arc = self.adjacency[here.index()][i];
                let bit = 1u64 << arc.edge;
                if used & bit != 0 { continue; }
                self.walk_edges.push(arc.edge);
                if arc.to == start && seen == all_nodes
                    && length.checked_sub(arc.overlap) == Some(self.target.len()) {
                    let proof = ReplayProof {
                        status: "PASS",
                        path_type: "circular",
                        used_segments: self.graph.names.len(),
                        used_links: self.walk_edges.len(),
                        states: self.states,
                        transitions: self.transitions,
                        starts: self.starts,
                    };
                    self.walk_edges.pop();
                    return Ok(Some(proof));
                }
                // The initial step is revisited only to close the circle. All L
                // edges remain legal, but each may be traversed at most once.
                if arc.to != start {
                    let next_seq = &self.graph.sequences[arc.to.node][usize::from(arc.to.reverse)];
                    let piece = &next_seq[arc.overlap..];
                    if let Some(end) = length.checked_add(piece.len()) {
                        if end <= self.target.len() + self.max_overlap
                            && position + end <= self.repeated_target.len()
                            && self.repeated_target[position + length..position + end] == *piece {
                            self.walk.push(arc.to);
                            let found = self.visit_existing(start, arc.to, position, end,
                                used | bit, seen | (1u64 << arc.to.node), all_nodes)?;
                            self.walk.pop();
                            if found.is_some() {
                                self.walk_edges.pop();
                                return Ok(found);
                            }
                        }
                    }
                }
                self.walk_edges.pop();
            }
            Ok(None)
        }

        fn prove_existing(mut self) -> Result<ReplayProof, OrgraftError> {
            let overlaps = self.graph.edges.iter().map(|edge| edge.old).collect::<Vec<_>>();
            if !self.p_valid(&overlaps) {
                return Err(invalid("an existing P path has no matching oriented L adjacency"));
            }
            if self.graph.sequences.iter().any(|pair| {
                pair[0].len() > self.target.len() + self.max_overlap
            }) {
                return Err(invalid("an S segment is longer than the accepted circle with maximum overlap"));
            }
            let anchor = (0..self.graph.names.len()).max_by_key(|&i| {
                (self.graph.sequences[i][0].len(), &self.graph.names[i])
            }).expect("at least one S");
            let all_nodes = if self.graph.names.len() == 64 {
                u64::MAX
            } else {
                (1u64 << self.graph.names.len()) - 1
            };
            for reverse in [false, true] {
                let start = Step { node: anchor, reverse };
                let query = &self.graph.sequences[anchor][usize::from(reverse)];
                let query_text = std::str::from_utf8(query).expect("validated ASCII DNA");
                let mut cursor = 0;
                while cursor < self.target.len() {
                    let offset = std::str::from_utf8(&self.repeated_target[cursor..])
                        .expect("validated ASCII DNA")
                        .find(query_text);
                    let Some(offset) = offset else { break };
                    let position = cursor + offset;
                    if position >= self.target.len() { break; }
                    self.starts += 1;
                    if self.starts > MAX_STARTS {
                        return Err(invalid("existing-circle start limit reached; human review required"));
                    }
                    self.walk.push(start);
                    let found = self.visit_existing(start, start, position, query.len(),
                        0, 1u64 << anchor, all_nodes)?;
                    self.walk.pop();
                    if let Some(proof) = found { return Ok(proof); }
                    cursor = position + 1;
                }
            }
            Err(invalid("accepted closed walk not found in the existing graph"))
        }

        fn record(&mut self) {
            let overlaps = self.assignments.iter().map(|n| n.expect("all L visited"))
                .collect::<Vec<_>>();
            if !self.p_valid(&overlaps) { return; }
            let identity = (overlaps.clone(), canonical_cycle(&self.walk, &self.walk_edges));
            if self.identities.contains(&identity) { return; }
            self.identities.push(identity);
            self.solutions.push(Solution {
                overlaps, walk: self.walk.clone(), walk_edges: self.walk_edges.clone(),
            });
        }

        fn visit(
            &mut self, start: Step, here: Step, position: usize, length: usize,
            used: u64, seen: u64, all_edges: u64, all_nodes: u64,
        ) -> Result<(), OrgraftError> {
            self.states += 1;
            if self.states > MAX_STATES {
                return Err(invalid("complete-walk search limit reached; human review required"));
            }
            for i in 0..self.adjacency[here.index()].len() {
                self.transitions += 1;
                if self.transitions > MAX_TRANSITIONS {
                    return Err(invalid("complete-walk transition limit reached; human review required"));
                }
                let arc = self.adjacency[here.index()][i];
                let bit = 1u64 << arc.edge;
                if used & bit != 0 { continue; }
                let next_used = used | bit;
                self.assignments[arc.edge] = Some(arc.overlap);
                self.walk_edges.push(arc.edge);
                if next_used == all_edges && arc.to == start
                    && seen == all_nodes && length.checked_sub(arc.overlap) == Some(self.target.len()) {
                    self.record();
                } else if next_used != all_edges {
                    let next_seq = &self.graph.sequences[arc.to.node][usize::from(arc.to.reverse)];
                    let piece = &next_seq[arc.overlap..];
                    if let Some(end) = length.checked_add(piece.len()) {
                        if end <= self.target.len() + self.max_overlap
                            && position + end <= self.repeated_target.len()
                            && self.repeated_target[position + length..position + end] == *piece {
                            self.walk.push(arc.to);
                            self.visit(start, arc.to, position, end, next_used,
                                seen | (1u64 << arc.to.node), all_edges, all_nodes)?;
                            self.walk.pop();
                        }
                    }
                }
                self.walk_edges.pop();
                self.assignments[arc.edge] = None;
                if self.solutions.len() > 1 { return Ok(()); }
            }
            Ok(())
        }

        fn prove(mut self) -> Result<Proof, OrgraftError> {
            let anchor = (0..self.graph.names.len()).max_by_key(|&i| {
                (self.graph.sequences[i][0].len(), &self.graph.names[i])
            }).expect("at least one S");
            let all_nodes = if self.graph.names.len() == 64 {
                u64::MAX
            } else {
                (1u64 << self.graph.names.len()) - 1
            };
            let all_edges = if self.graph.edges.len() == 64 {
                u64::MAX
            } else {
                (1u64 << self.graph.edges.len()) - 1
            };
            for reverse in [false, true] {
                let start = Step { node: anchor, reverse };
                let query = &self.graph.sequences[anchor][usize::from(reverse)];
                let query_text = std::str::from_utf8(query).expect("validated ASCII DNA");
                let mut cursor = 0;
                while cursor < self.target.len() {
                    let offset = std::str::from_utf8(&self.repeated_target[cursor..])
                        .expect("validated ASCII DNA")
                        .find(query_text);
                    let Some(offset) = offset else { break };
                    let position = cursor + offset;
                    if position >= self.target.len() { break; }
                    self.starts += 1;
                    if self.starts > MAX_STARTS {
                        return Err(invalid("accepted-circle start limit reached; human review required"));
                    }
                    self.walk.push(start);
                    self.visit(start, start, position, query.len(), 0, 1u64 << anchor,
                        all_edges, all_nodes)?;
                    self.walk.pop();
                    if self.solutions.len() > 1 {
                        return Err(invalid("multiple valid circular walks or overlap vectors; human review required"));
                    }
                    cursor = position + 1;
                }
            }
            if self.solutions.len() != 1 {
                return Err(invalid("no complete accepted circular walk uses every S and L exactly once"));
            }
            Ok(Proof {
                solution: self.solutions.pop().expect("one solution"),
                states: self.states,
                transitions: self.transitions,
                starts: self.starts,
                candidate_counts: self.graph.edges.iter().map(|edge| edge.candidates.len()).collect(),
            })
        }
    }

    pub(super) fn validate_existing_circle(
        candidate_gfa: &[u8], accepted: &str,
    ) -> Result<ReplayProof, OrgraftError> {
        if candidate_gfa.len() > MAX_GFA_BYTES {
            return Err(invalid("existing graph exceeds bounded replay input size"));
        }
        if accepted.is_empty() || accepted.len() > MAX_TARGET_BP
            || !accepted.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b)) {
            return Err(invalid("one nonempty IUPAC-DNA accepted sequence within the size limit is required"));
        }
        let text = std::str::from_utf8(candidate_gfa)
            .map_err(|_| invalid("existing graph is not UTF-8"))?;
        let graph = parse_graph(text, false)?;
        Search::new(&graph, accepted.as_bytes()).prove_existing()
    }

    pub(super) fn validate_linear_singleton(
        graph: &Gfa, accepted: &str,
    ) -> Result<ReplayProof, OrgraftError> {
        if !graph.links.is_empty() || graph.segments.len() != 1 || graph.order.len() != 1 {
            return Err(invalid("a link-free accepted linear model must contain exactly one S segment"));
        }
        if accepted.is_empty() || accepted.len() > MAX_TARGET_BP
            || !accepted.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b)) {
            return Err(invalid("one nonempty IUPAC-DNA accepted sequence within the size limit is required"));
        }
        let node = &graph.order[0];
        let segment = &graph.segments[node];
        let sequence = segment.sequence.to_ascii_uppercase();
        if !sequence.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b)) {
            return Err(invalid("linear S node contains non-IUPAC DNA"));
        }
        if sequence.len() != accepted.len() {
            return Err(invalid("the existing linear S node and accepted FASTA differ in length"));
        }
        let doubled = format!("{accepted}{accepted}");
        if !doubled.contains(&sequence)
            && !doubled.contains(&reverse_complement_iupac(&sequence)) {
            return Err(invalid("LINEAR_SINGLETON_SEQUENCE_DISAGREEMENT: the existing linear S node is not circular-frame equivalent to the accepted FASTA"));
        }
        for row in &graph.other_lines {
            if row.first().map(String::as_str) != Some("P") || row.len() < 4
                || row[3] != "*" || row[2].split(',').count() != 1 {
                return Err(invalid("link-free linear graph has an unsupported or illegal path record"));
            }
            let token = &row[2];
            if token != &format!("{node}+") && token != &format!("{node}-") {
                return Err(invalid("link-free linear P path does not name its only S node"));
            }
        }
        Ok(ReplayProof {
            status: "REVIEW_linear_no_closure",
            path_type: "linear_singleton_sequence_equivalent_no_closing_L",
            used_segments: 1,
            used_links: 0,
            states: 0,
            transitions: 0,
            starts: 1,
        })
    }

    fn patch_lines(text: &str, graph: &Graph, overlaps: &[usize]) -> Result<String, OrgraftError> {
        let mut by_line = BTreeMap::new();
        for (edge, overlap) in graph.edges.iter().zip(overlaps) {
            by_line.insert(edge.line_index, *overlap);
        }
        let mut output = String::with_capacity(text.len());
        for (line_index, line) in text.split_inclusive('\n').enumerate() {
            let Some(&overlap) = by_line.get(&line_index) else {
                output.push_str(line);
                continue;
            };
            let (content, ending) = if let Some(content) = line.strip_suffix("\r\n") {
                (content, "\r\n")
            } else if let Some(content) = line.strip_suffix('\n') {
                (content, "\n")
            } else {
                (line, "")
            };
            let mut fields = content.split('\t').collect::<Vec<_>>();
            if fields.len() < 6 || fields[0] != "L" {
                return Err(invalid("source L line changed during overlap repair"));
            }
            let new_overlap = format!("{overlap}M");
            fields[5] = &new_overlap;
            output.push_str(&fields.join("\t"));
            output.push_str(ending);
        }
        Ok(output)
    }

    fn unchanged_except_l_overlap(source: &str, patched: &str) -> bool {
        let before = source.split_inclusive('\n').collect::<Vec<_>>();
        let after = patched.split_inclusive('\n').collect::<Vec<_>>();
        if before.len() != after.len() { return false; }
        before.iter().zip(after.iter()).all(|(old, new)| {
            if !body(old).starts_with("L\t") { return old == new; }
            let old_fields = body(old).split('\t').collect::<Vec<_>>();
            let new_fields = body(new).split('\t').collect::<Vec<_>>();
            old_fields.len() == new_fields.len()
                && old_fields[..5] == new_fields[..5]
                && old_fields[6..] == new_fields[6..]
                && old.ends_with("\r\n") == new.ends_with("\r\n")
                && old.ends_with('\n') == new.ends_with('\n')
        })
    }

    fn repair(text: &str, accepted: &str) -> Result<(String, Proof), OrgraftError> {
        if accepted.is_empty() || accepted.len() > MAX_TARGET_BP
            || !accepted.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b)) {
            return Err(invalid("one nonempty IUPAC-DNA accepted sequence within the size limit is required"));
        }
        let graph = parse_graph(text, true)?;
        let proof = Search::new(&graph, accepted.as_bytes()).prove()?;
        let patched = patch_lines(text, &graph, &proof.solution.overlaps)?;
        if !unchanged_except_l_overlap(text, &patched) {
            return Err(invalid("overlap repair changed a protected GFA field"));
        }
        // Parse again so any accidental non-L mutation or invalid nM cannot escape.
        let checked = parse_graph(&patched, true)?;
        if checked.names != graph.names || checked.paths.len() != graph.paths.len()
            || checked.edges.len() != graph.edges.len() {
            return Err(invalid("overlap repair changed the graph contract"));
        }
        Ok((patched, proof))
    }

    pub(super) fn run(
        gfa_path: &Path, fasta_path: &Path, out_dir: &Path, subgraph: &str,
    ) -> Result<(), OrgraftError> {
        if out_dir.exists() {
            return Err(invalid("output directory already exists; overwrite is forbidden"));
        }
        let source_gfa = fs::read(gfa_path)?;
        if source_gfa.len() > MAX_GFA_BYTES {
            return Err(invalid("final GFA exceeds the bounded overlap-only input size"));
        }
        let source_fasta = fs::read(fasta_path)?;
        let text = std::str::from_utf8(&source_gfa)
            .map_err(|_| invalid("final GFA must be UTF-8"))?;
        let records = read_fasta(fasta_path)?;
        if records.len() != 1 {
            return Err(invalid("final accepted FASTA must contain exactly one sequence"));
        }
        let (patched, proof) = repair(text, &records[0].sequence)?;
        if fs::read(gfa_path)? != source_gfa || fs::read(fasta_path)? != source_fasta {
            return Err(invalid("input changed while proving the overlap-only repair"));
        }
        fs::create_dir(out_dir)?;
        let graph_dir = out_dir.join(subgraph);
        let logs_dir = out_dir.join("logs");
        fs::create_dir(&graph_dir)?;
        fs::create_dir(&logs_dir)?;
        let graph_out = graph_dir.join(format!("rebuild_{subgraph}.gfa"));
        let fasta_out = graph_dir.join(format!("rebuild_{subgraph}.fasta"));
        fs::write(&graph_out, patched.as_bytes())?;
        fs::write(&fasta_out, &source_fasta)?;
        let old = parse_graph(text, true)?;
        let changed = old.edges.iter().zip(&proof.solution.overlaps)
            .filter(|(edge, overlap)| edge.old != **overlap).count();
        let mut report = format!(
            "section\tkey\tvalue\nrun\tstatus\tPASS_unique_complete_circle\ninput\tfinal_gfa_sha256\t{}\ninput\taccepted_fasta_sha256\t{}\noutput\tfinal_gfa_sha256\t{}\noutput\taccepted_fasta_sha256\t{}\nproof\tused_segments\t{}\nproof\tused_links\t{}\nproof\tP_records_preserved\t{}\nproof\tread_support\tnot_recomputed\nproof\tsearch_states\t{}\nproof\tsearch_starts\t{}\nproof\tchanged_L_overlaps\t{}\nproof\tselected_walk\t{}\nproof\tselected_L_order\t{}\nproof\tselected_overlaps\t{}\nproof\tcandidate_counts\t{}\n",
            sha256_hex(&source_gfa), sha256_hex(&source_fasta),
            sha256_hex(patched.as_bytes()), sha256_hex(&source_fasta),
            old.names.len(), proof.solution.walk_edges.len(), old.paths.len(),
            proof.states, proof.starts, changed,
            proof.solution.walk.iter().map(|step| format!("{}{}", old.names[step.node], if step.reverse { '-' } else { '+' }))
                .collect::<Vec<_>>().join(","),
            proof.solution.walk_edges.iter().map(|index| (index + 1).to_string()).collect::<Vec<_>>().join(","),
            proof.solution.overlaps.iter().map(|n| format!("{n}M")).collect::<Vec<_>>().join(","),
            proof.candidate_counts.iter().map(usize::to_string).collect::<Vec<_>>().join(","),
        );
        report.push_str(&format!("proof\tsearch_transitions\t{}\n", proof.transitions));
        fs::write(logs_dir.join(format!("rebuild_{subgraph}_overlap_only.tsv")), report)?;
        Ok(())
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::commands::rebuild::RebuildOptions;

        const TWO_NODE: &str = "H\tVN:Z:1.0\nS\ta\tAAC\tSN:Z:first\nS\tb\tCGT\nL\ta\t+\tb\t+\t0M\tXX:Z:keep\nL\tb\t+\ta\t+\t0M\nP\tp\ta+,b+\t*\tPP:Z:keep\n";

        #[test]
        fn repairs_only_l_overlap_and_keeps_every_other_field() {
            let (patched, proof) = repair(TWO_NODE, "AACGT").unwrap();
            assert_eq!(proof.solution.overlaps, vec![1, 0]);
            assert_eq!(proof.candidate_counts, vec![2, 1]);
            assert_eq!(patched, TWO_NODE.replacen(
                "L\ta\t+\tb\t+\t0M\tXX:Z:keep",
                "L\ta\t+\tb\t+\t1M\tXX:Z:keep", 1,
            ));
        }

        #[test]
        fn closing_overlap_of_a_single_node_is_proven() {
            let source = "S\ta\tACGTAC\nL\ta\t+\ta\t+\t0M\n";
            let (patched, proof) = repair(source, "ACGT").unwrap();
            assert_eq!(proof.solution.overlaps, vec![2]);
            assert!(patched.contains("L\ta\t+\ta\t+\t2M"));
        }

        #[test]
        fn a_complete_walk_may_revisit_a_repeat_node() {
            let source = "S\ta\tAAC\nS\tb\tCTT\nL\ta\t+\ta\t+\t0M\nL\ta\t+\tb\t+\t0M\nL\tb\t+\ta\t+\t0M\n";
            let (_, proof) = repair(source, "AACAACTT").unwrap();
            assert_eq!(proof.solution.overlaps, vec![0, 1, 0]);
            assert_eq!(proof.solution.walk_edges.len(), 3);
            assert_eq!(proof.solution.walk.len(), 3);
        }

        #[test]
        fn reverse_oriented_links_are_spelled_in_the_correct_direction() {
            let source = "S\ta\tAAC\nS\tb\tTCG\nL\ta\t+\tb\t-\t0M\nL\tb\t-\ta\t+\t0M\n";
            let (_, proof) = repair(source, "AACG").unwrap();
            assert_eq!(proof.solution.overlaps, vec![1, 1]);
        }

        #[test]
        fn distinct_valid_overlap_vectors_fail_closed() {
            let source = "S\ta\tAAAA\nS\tb\tAAAA\nL\ta\t+\tb\t+\t0M\nL\tb\t+\ta\t+\t0M\n";
            let error = repair(source, "AAAAAAA").err().unwrap().to_string();
            assert!(error.contains("multiple valid"), "{error}");
        }

        #[test]
        fn disconnected_final_graph_cannot_be_inferred() {
            let source = "S\ta\tAAC\nS\tb\tCGT\nL\ta\t+\tb\t+\t0M\n";
            let error = repair(source, "AACGT").err().unwrap().to_string();
            assert!(error.contains("no complete accepted"), "{error}");
        }

        #[test]
        fn explicit_p_overlap_is_a_constraint() {
            let source = TWO_NODE.replace("P\tp\ta+,b+\t*", "P\tp\ta+,b+\t0M");
            let error = repair(&source, "AACGT").err().unwrap().to_string();
            assert!(error.contains("no complete accepted"), "{error}");
        }

        #[test]
        fn unsupported_path_record_fails_closed() {
            let source = format!("{TWO_NODE}W\twalk\t0\ta\t0\t3\t>\n");
            let error = repair(&source, "AACGT").err().unwrap().to_string();
            assert!(error.contains("unsupported GFA record"), "{error}");
        }

        #[test]
        fn cli_requires_fresh_explicit_output_and_rejects_force() {
            let base = [
                "--edited-gfa", "graph.gfa", "--polished-fasta", "accepted.fasta",
                "--overlap-only",
            ];
            let without_out = base.iter().map(|s| s.to_string()).collect::<Vec<_>>();
            assert!(RebuildOptions::from_args(&without_out).is_err());
            let with_force = base.iter().chain(["--out-dir", "new-output", "--force"].iter())
                .map(|s| s.to_string()).collect::<Vec<_>>();
            assert!(RebuildOptions::from_args(&with_force).is_err());
            let valid = base.iter().chain(["--out-dir", "new-output"].iter())
                .map(|s| s.to_string()).collect::<Vec<_>>();
            assert!(RebuildOptions::from_args(&valid).unwrap().overlap_only);
        }

        #[test]
        fn fresh_output_copies_fasta_bytes_and_never_overwrites() {
            let stamp = std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
            let base = std::env::temp_dir().join(format!(
                "orgraft-overlap-only-{}-{stamp}", std::process::id(),
            ));
            fs::create_dir(&base).unwrap();
            let gfa = base.join("input.gfa");
            let fasta = base.join("accepted.fasta");
            let out = base.join("output");
            fs::write(&gfa, TWO_NODE).unwrap();
            let source_fasta = b">accepted carefully annotated\nAACGT\n";
            fs::write(&fasta, source_fasta).unwrap();
            run(&gfa, &fasta, &out, "subgraph_001").unwrap();
            assert_eq!(fs::read(out.join("subgraph_001/rebuild_subgraph_001.fasta")).unwrap(), source_fasta);
            assert_eq!(fs::read_to_string(out.join("subgraph_001/rebuild_subgraph_001.gfa")).unwrap(),
                TWO_NODE.replacen("L\ta\t+\tb\t+\t0M", "L\ta\t+\tb\t+\t1M", 1));
            assert!(run(&gfa, &fasta, &out, "subgraph_001").is_err());
            assert_eq!(fs::read(&gfa).unwrap(), TWO_NODE.as_bytes());
            let _ = fs::remove_dir_all(&base);
        }

        #[test]
        fn existing_circle_accepts_a_repeat_alternative_without_using_every_link() {
            let source = "S\ta\tAAC\nS\tb\tCGT\nS\tc\tTAA\nL\ta\t+\tb\t+\t1M\nL\tb\t+\tc\t+\t1M\nL\tc\t+\ta\t+\t1M\nL\ta\t+\tc\t+\t0M\nP\tp\ta+,b+,c+\t*\n";
            let proof = validate_existing_circle(source.as_bytes(), "AACGTA").unwrap();
            assert_eq!(proof.used_segments, 3);
            assert_eq!(proof.used_links, 3);
        }

        #[test]
        fn existing_circle_rejects_wrong_direction_even_if_nodes_tile() {
            let source = "S\ta\tAAC\nS\tb\tCGT\nS\tc\tTAA\nL\ta\t+\tc\t+\t0M\nL\tc\t+\tb\t+\t0M\nL\tb\t+\ta\t+\t0M\n";
            let error = validate_existing_circle(source.as_bytes(), "AACGTA")
                .err().unwrap().to_string();
            assert!(error.contains("accepted closed walk not found"), "{error}");
        }

        #[test]
        fn existing_circle_rejects_illegal_unused_link_or_p_path() {
            let source = "S\ta\tAAC\nS\tb\tCGT\nL\ta\t+\tb\t+\t1M\nL\tb\t+\ta\t+\t0M\n";
            assert!(validate_existing_circle(source.as_bytes(), "AACGT").is_ok());
            let bad_link = format!("{source}L\ta\t+\tb\t+\t2M\n");
            assert!(validate_existing_circle(bad_link.as_bytes(), "AACGT")
                .err().unwrap().to_string().contains("existing L overlap"));
            let bad_path = format!("{source}P\tp\tb+,b+\t*\n");
            assert!(validate_existing_circle(bad_path.as_bytes(), "AACGT")
                .err().unwrap().to_string().contains("P path"));
        }

        #[test]
        fn link_free_linear_singleton_replays_without_requiring_a_circle() {
            let graph = super::super::Gfa {
                headers: Vec::new(),
                segments: [("a".to_string(), super::super::Segment {
                    name: "a".to_string(), sequence: "AACGT".to_string(), tags: Vec::new(),
                })].into(),
                order: vec!["a".to_string()],
                links: Vec::new(),
                other_lines: vec![vec!["P".to_string(), "p".to_string(),
                    "a+".to_string(), "*".to_string()]],
            };
            assert_eq!(validate_linear_singleton(&graph, "AACGT").unwrap().status,
                "REVIEW_linear_no_closure");
            assert_eq!(validate_linear_singleton(&graph, "CGTAA").unwrap().status,
                "REVIEW_linear_no_closure");
            assert!(validate_linear_singleton(&graph, "AACGA")
                .err().unwrap().to_string().contains("not circular-frame equivalent"));
            let mut illegal = graph.clone();
            illegal.other_lines[0][2] = "a+,a+".to_string();
            assert!(validate_linear_singleton(&illegal, "AACGT").is_err());
        }
    }
}


const HELP: &str = r#"orgraft rebuild

Rebuild a compact verified GFA from an edited GFA and polished FASTA.

Usage:
  orgraft rebuild --organelle NAME --subgraph ID --edited-gfa FILE --polished-fasta FILE [options]

Inputs:
  --organelle NAME              organelle name for this rebuild run [mito]
  --subgraph ID                 subgraph/ring id [subgraph_001]
  --edited-gfa FILE             edited GFA; original resolve component with --resolve-details
  --polished-fasta FILE         final polished FASTA
  --checked-gfa FILE            original checked graph; requires --resolve-details
  --resolve-details FILE        recorded resolve walk; enables topology and path repair
  --ordinary-source-dir DIR     validated ordinary round; replay its source transition only if the recorded path fails
  --joint-ir-trial-gfa FILE      new joint IR candidate GFA; requires paired evidence
  --joint-ir-evidence FILE       hash-bound joint IR table from postrepair
  --overlap-only                repair only L overlaps in a final accepted circular GFA
  --structure-corrections FILE  hash-bound, evidence-backed graph corrections (optional)
  --reads FILE                  organelle reads for repeat-path support
  --soft-paths FILE             tool paths file [soft_paths.txt]

Outputs:
  --out-dir DIR                 rebuild root output directory [results/rebuild]
  --force                       replace only this component's outputs

Additional Parameters:
  --threads N                   minimap2 threads for node-to-polished projection [4]
  --image-reference-fasta FILE  reference FASTA for graph colouring; enables PDF/SVG export
  --repeat-resolution-node NODE constrain one repeat using complete spanning-read support

Layout: OUT/SUBGRAPH/rebuild_SUBGRAPH* plus OUT/logs/*.tsv (pdf/svg need reference)
"#;

const DEFAULT_SOFT_PATHS: &str = "soft_paths.txt";
const DEFAULT_THREADS: usize = 4;
const DEFAULT_ORGANELLE: &str = "mito";
const DEFAULT_SUBGRAPH: &str = "subgraph_001";
const DEFAULT_OUT_DIR: &str = "results/rebuild";
const REPEAT_PATH_TEMPLATE_FLANK_BP: usize = 2_000;
const REPEAT_PATH_REQUIRED_FLANK_BP: usize = 1_000;
const REPEAT_PATH_ALIGNMENT_TOLERANCE_BP: usize = 50;
const REPEAT_PATH_MIN_IDENTITY: f64 = 0.95;
const REPEAT_PATH_SECONDARY_LIMIT: usize = 50;
// A topology-declared repeat may be shorter than 1 kb. Self-BLAST is only a
// candidate generator; two independent node placements must still anchor it.
const REPEAT_CORE_MIN_ALIGNMENT_BP: usize = 500;
const REPEAT_CORE_DEFAULT_ALIGNMENT_BP: usize = 1_000;
const REPEAT_ASSIGNMENT_MIN_ALIGNED_FRACTION: f64 = 0.75;
const REPEAT_ASSIGNMENT_MIN_SUBSTITUTION_IDENTITY: f64 = 0.99;

pub fn run(args: &[String]) -> Result<(), OrgraftError> {
    if args.is_empty() {
        println!("{HELP}");
        return Ok(());
    }
    if args.iter().any(|arg| arg == "-h" || arg == "--help") {
        println!("{HELP}");
        return Ok(());
    }
    if args.iter().any(|arg| arg == "--contract") {
        print_contract(&contract());
        return Ok(());
    }
    let options = RebuildOptions::from_args(args)?;
    if options.overlap_only {
        return overlap_only::run(
            &options.edited_gfa,
            &options.polished_fasta,
            &options.out_dir,
            &options.subgraph,
        );
    }
    run_rebuild(&options)
}

fn joint_ir_trial_sequences(
    trial_path: &Path,
    evidence_path: &Path,
    checked_path: &Path,
    details_path: &Path,
    polished_path: &Path,
    contract: &ResolveContract,
) -> Result<(BTreeMap<String, String>, bool), OrgraftError> {
    let fail = |message: &str| OrgraftError::InvalidArgument(format!("joint IR evidence: {message}"));
    let text = fs::read_to_string(evidence_path)?;
    if text.len() > 100_000 {
        return Err(fail("event table is too large"));
    }
    let mut rows = BTreeMap::new();
    for (index, line) in text.lines().enumerate() {
        if index == 0 {
            if line != "key\tvalue" {
                return Err(fail("event table header changed"));
            }
            continue;
        }
        let (key, value) = line.split_once('\t').ok_or_else(|| fail("malformed event row"))?;
        if rows.insert(key.to_string(), value.to_string()).is_some() {
            return Err(fail("duplicate event key"));
        }
    }
    let get = |key: &str| -> Result<&str, OrgraftError> {
        rows.get(key).map(String::as_str).ok_or_else(|| fail("required event key is absent"))
    };
    for (key, path) in [
        ("checked_gfa_sha256", checked_path),
        ("resolve_details_sha256", details_path),
        ("candidate_fasta_sha256", polished_path),
        ("candidate_gfa_sha256", trial_path),
    ] {
        if get(key)? != file_sha(path)? {
            return Err(fail("input hash differs from the joint IR event table"));
        }
    }
    if get("status")? != "selected_for_ordinary_validation"
        || get("R2_path_links")? != "verified_zero_overlap"
        || get("opposite_pair_count")? != "1"
    {
        return Err(fail("joint event was not selected with verified source links"));
    }
    let count = |key: &str| -> Result<usize, OrgraftError> {
        get(key)?.parse().map_err(|_| fail("invalid read count"))
    };
    let pool = count("fixed_read_pool")?;
    let main = count("main_ALT")?;
    let other = count("other_complete_ALT")?;
    let remainder = count("remainder")?;
    if main <= remainder
        || main.checked_add(other).and_then(|v| v.checked_add(remainder)) != Some(pool)
        || get("strict_formula")? != format!("{main} > {remainder}")
    {
        return Err(fail("fixed read pool fails the unchanged strict gate"));
    }
    let repeat = get("repeat_node")?;
    if contract.walk.iter().filter(|step| step.source_node == repeat).count() != 2 {
        return Err(fail("joint repeat is not present twice in the recorded walk"));
    }
    let trial = Gfa::read(trial_path)?;
    if trial.segments.keys().cloned().collect::<BTreeSet<_>>()
        != contract.source_sequences.keys().cloned().collect()
    {
        return Err(fail("trial graph changed source node identities"));
    }
    let trial_links = trial.links.iter().map(|link| {
        (transition_key(&link.from_name, link.from_orient, &link.to_name, link.to_orient), link.overlap.clone())
    }).collect::<BTreeSet<_>>();
    let source_links = contract.topology_links.iter().map(|link| {
        (physical_key(link), link.overlap.clone())
    }).collect::<BTreeSet<_>>();
    if trial.links.len() != trial_links.len()
        || contract.topology_links.len() != source_links.len()
        || trial_links != source_links
    {
        return Err(fail("trial graph changed the checked connections or overlaps"));
    }
    let mut sequences = trial.segments.iter().map(|(name, node)| (name.clone(), node.sequence.clone())).collect::<BTreeMap<_, _>>();
    if sequences[repeat].len() <= contract.source_sequences[repeat].len() {
        return Err(fail("selected shared IR node did not grow"));
    }
    // Joint polish may emit the two single-copy arms in the reference FASTA
    // direction while the recorded resolve walk uses their opposite strand.
    // Determine the arm basis from two internal checked-node anchors.
    // The shared IR node retains its declared orientation in either case.
    let mut arm_bases = Vec::new();
    for (name, source) in &contract.source_sequences {
        if name == repeat {
            continue;
        }
        if source.len() < 2200 {
            return Err(fail("single-copy arm is too short for orientation anchors"));
        }
        let candidate = &sequences[name];
        let opposite = reverse_complement(candidate);
        let seeds = [&source[1000..1100], &source[source.len() - 1100..source.len() - 1000]];
        let direct = seeds.iter().all(|seed| candidate.contains(seed));
        let reverse = seeds.iter().all(|seed| opposite.contains(seed));
        if direct == reverse {
            return Err(fail("single-copy arm orientation is ambiguous or lacks checked anchors"));
        }
        arm_bases.push(reverse);
    }
    if arm_bases.is_empty() || arm_bases.iter().any(|basis| *basis != arm_bases[0]) {
        return Err(fail("single-copy arms disagree on orientation basis"));
    }
    let normalized = arm_bases[0];
    if normalized {
        for (name, sequence) in &mut sequences {
            if name != repeat {
                *sequence = reverse_complement(sequence);
            }
        }
    }
    Ok((sequences, normalized))
}

fn recorded_candidate_plan(graph: &mut Gfa, target: &FastaRecord, fasta: &Path,
    contract: &ResolveContract, minimap2: &Path, blastn: &Path, options: &RebuildOptions,
    paths: &OutputPaths) -> Result<BuildResult,OrgraftError> {
    let original_graph = graph.clone();
    Ok(match build_plan(
            graph,
            target,
            contract,
            minimap2,
            options.threads,
            paths,
        ) {
            Ok(result) => result,
            Err(path_error) => {
                // Both generators use the original source graph. Existing rebuilt
                // models are comparison evidence and are never reconstruction inputs.
                *graph = original_graph;
                let alternative = (|| {
                    let mapping = run_node_projection(
                        &*graph,
                        target,
                        minimap2,
                        options.threads,
                        paths,
                        options.keep_debug,
                    )?;
                    let classes = classify_nodes(&*graph, &mapping);
                    let proposal = infer_verified_sequences_from_repeat_cores(
                        &*graph,
                        target,
                        fasta,
                        &mapping,
                        &classes,
                        blastn,
                        paths,
                        options.keep_debug,
                    )?;
                    accept_node_sequences(
                        graph,
                        target,
                        contract,
                        proposal.sequences,
                        paths,
                        "node_projection_with_exact_graph_replay",
                    )
                })();
                match alternative {
                    Ok((plan, classes, mut summary)) => {
                        summary.push(("whole_path_generator_diagnostic", path_error.to_string().chars().take(2000).collect()));
                        (plan, classes, summary)
                    }
                    Err(projection_error) => return Err(OrgraftError::InvalidArgument(format!(
                        "{}{}",
                        reconstruction_failure_tag(&projection_error),
                        format_args!("no verified reconstruction; whole-path generator: {path_error}; node-projection generator: {projection_error}")
                    ))),
                }
            }
        })
}

fn recorded_candidate_or_source(graph: &mut Gfa, target: &FastaRecord, fasta: &Path,
    contract: &mut ResolveContract, minimap2: &Path, blastn: &Path, options: &RebuildOptions,
    paths: &OutputPaths) -> Result<BuildResult,OrgraftError> {
    let original_graph=graph.clone();
    match recorded_candidate_plan(graph,target,fasta,contract,minimap2,blastn,options,paths) {
        Ok(result)=>Ok(result),
        Err(original_error) if options.ordinary_source_dir.is_some()=>{
            // Consult ordinary history only after both original generators fail.
            let derived=ordinary_source::produce(options,paths)?;
            *graph=original_graph;
            let mut replay_paths=paths.clone();
            replay_paths.debug_dir=paths.debug_dir.join("ordinary_replay");
            let result=recorded_candidate_plan(graph,target,fasta,&derived,minimap2,blastn,options,&replay_paths)
                .map_err(|e| OrgraftError::InvalidArgument(format!("ordinary source transition verified, but final reconstruction failed: {e}; original path: {original_error}")))?;
            *contract=derived;Ok(result)
        }
        Err(error)=>Err(error),
    }
}

fn run_rebuild(options: &RebuildOptions) -> Result<(), OrgraftError> {
    let started = Instant::now();
    let edited_gfa = canonicalize_existing(&options.edited_gfa, "--edited-gfa")?;
    let polished_fasta = canonicalize_existing(&options.polished_fasta, "--polished-fasta")?;
    let reads = options
        .reads
        .as_ref()
        .map(|path| canonicalize_existing(path, "--reads"))
        .transpose()?;
    if options.repeat_resolution_node.is_some() && reads.is_none() {
        return Err(OrgraftError::InvalidArgument(
            "--repeat-resolution-node requires --reads for complete flank-repeat-flank support"
                .to_string(),
        ));
    }
    let image_reference_fasta = options
        .image_reference_fasta
        .as_ref()
        .map(|path| canonicalize_existing(path, "--image-reference-fasta"))
        .transpose()?;
    let merged_template = options
        .merged_gfa_template
        .as_ref()
        .map(|path| canonicalize_existing(path, "--merged-gfa-template"))
        .transpose()?;
    let mut recorded_contract = match (&options.checked_gfa, &options.resolve_details) {
        (Some(checked), Some(details)) => Some(load_contract(
            &canonicalize_existing(checked, "--checked-gfa")?,
            &edited_gfa,
            &canonicalize_existing(details, "--resolve-details")?,
            &options.subgraph,
        )?),
        _ => None,
    };
    let joint_ir_sequences = if let (Some(trial), Some(evidence), Some(contract)) = (
        &options.joint_ir_trial_gfa,
        &options.joint_ir_evidence,
        &recorded_contract,
    ) {
        Some(joint_ir_trial_sequences(
            &canonicalize_existing(trial, "--joint-ir-trial-gfa")?,
            &canonicalize_existing(evidence, "--joint-ir-evidence")?,
            &canonicalize_existing(options.checked_gfa.as_ref().unwrap(), "--checked-gfa")?,
            &canonicalize_existing(options.resolve_details.as_ref().unwrap(), "--resolve-details")?,
            &polished_fasta,
            contract,
        )?)
    } else {
        None
    };
    // Manual structure corrections are TEMPORARILY DISABLED (decision 2026-09-12): the
    // CLI rejects --structure-corrections, so this binding is always None and the block
    // below is dormant. The full implementation is retained in rebuild_structure.rs and
    // can be re-enabled by restoring CLI parsing above. Detection of a source-graph
    // structural error (scenario C) happens in the both-generators-failed branch instead.
    let structural_corrections = if let Some(record) = &options.structure_corrections {
        let checked =
            canonicalize_existing(options.checked_gfa.as_ref().unwrap(), "--checked-gfa")?;
        let details = canonicalize_existing(
            options.resolve_details.as_ref().unwrap(),
            "--resolve-details",
        )?;
        Some(apply(
            &canonicalize_existing(record, "--structure-corrections")?,
            recorded_contract.as_mut().unwrap(),
            &[
                ("checked_gfa", &checked),
                ("edited_gfa", &edited_gfa),
                ("resolve_details", &details),
                ("polished_fasta", &polished_fasta),
            ],
        )?)
    } else {
        None
    };
    if let Some(contract) = &recorded_contract {
        if !contract.path_issues.is_empty() {
            return Err(OrgraftError::InvalidArgument(format!(
                "recorded resolve walk is not a legal graph path: {}; no rebuild written",
                contract.path_issues.join("; ")
            )));
        }
    }

    let paths = OutputPaths::new(&options.out_dir, &options.subgraph);
    let prior_repeat_resolution = paths.prepare(options)?;

    let raw_gfa = Gfa::read(&edited_gfa)?;
    let polished_records = read_fasta(&polished_fasta)?;
    if polished_records.len() != 1 {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} must contain exactly one FASTA record, found {}",
            polished_fasta.display(),
            polished_records.len()
        )));
    }
    let polished_record = &polished_records[0];
    let (mut merged_gfa, merge_mode) = if let Some(contract) = &recorded_contract {
        let mut graph = raw_gfa.clone();
        graph
            .order
            .retain(|node| contract.source_sequences.contains_key(node));
        graph
            .segments
            .retain(|node, _| contract.source_sequences.contains_key(node));
        for (node, sequence) in &contract.source_sequences {
            // A split_node correction removes the parent and introduces child nodes
            // that are absent from the raw merged GFA: create their segments/order
            // entries instead of assuming every corrected node already exists.
            if let Some(segment) = graph.segments.get_mut(node) {
                segment.sequence = sequence.clone();
            } else {
                graph.segments.insert(
                    node.clone(),
                    Segment {
                        name: node.clone(),
                        sequence: sequence.clone(),
                        tags: Vec::new(),
                    },
                );
                graph.order.push(node.clone());
            }
        }
        if structural_corrections.is_some() {
            if graph.other_lines.iter().any(|row| {
                matches!(
                    row.first().map(String::as_str),
                    Some("C" | "J" | "F" | "E" | "G" | "O" | "U")
                )
            }) {
                return Err(OrgraftError::InvalidArgument("structure corrections do not support additional unparsed GFA connection records".to_string()));
            }
            graph
                .other_lines
                .retain(|row| !matches!(row.first().map(String::as_str), Some("P" | "W")));
        }
        graph.links = contract
            .topology_links
            .iter()
            .map(|edge| Link {
                from_name: edge.from.clone(),
                from_orient: edge.from_orientation,
                to_name: edge.to.clone(),
                to_orient: edge.to_orientation,
                overlap: edge.overlap.clone(),
                tags: edge.tags.clone(),
            })
            .collect();
        (
            graph,
            if structural_corrections.is_some() {
                "checked_graph_with_evidence_backed_structure_corrections"
            } else {
                "checked_compaction_and_recorded_resolve_path"
            }
            .to_string(),
        )
    } else if let Some(template) = &merged_template {
        (Gfa::read(template)?, "template".to_string())
    } else {
        merge_unambiguous_gfa(&raw_gfa)
    };

    if options.keep_debug {
        merged_gfa.write(&paths.debug_dir.join("merged_skeleton.gfa"))?;
    }

    let soft_paths = read_soft_paths(&options.soft_paths)?;
    let minimap2 = options
        .minimap2
        .clone()
        .or_else(|| soft_paths.get("minimap2").cloned())
        .unwrap_or_else(|| PathBuf::from("minimap2"));
    let blastn = options
        .blastn
        .clone()
        .or_else(|| soft_paths.get("blastn").cloned())
        .unwrap_or_else(|| PathBuf::from("blastn"));
    let support_source_gfa = merged_gfa.clone();
    let (sequence_plan, node_classes, mut path_summary) = if let (Some(contract), Some((sequences, normalized))) =
        (&recorded_contract, joint_ir_sequences)
    {
        let (plan, classes, mut summary) = accept_node_sequences(
            &mut merged_gfa,
            polished_record,
            contract,
            sequences,
            &paths,
            "joint_ir_hash_bound_original_checked_path",
        )?;
        summary.push(("joint_ir_evidence", "PASS_hash_bound_positive_gate".to_string()));
        summary.push(("joint_ir_arms_reverse_normalized", normalized.to_string()));
        (plan, classes, summary)
    } else if let Some(contract) = recorded_contract.as_mut() {
        recorded_candidate_or_source(&mut merged_gfa,polished_record,&polished_fasta,contract,&minimap2,&blastn,options,&paths)?
    } else {
        let mapping = run_node_projection(
            &merged_gfa,
            polished_record,
            &minimap2,
            options.threads,
            &paths,
            options.keep_debug,
        )?;

        let node_classes = classify_nodes(&merged_gfa, &mapping);
        let sequence_plan = infer_verified_sequences_from_repeat_cores(
            &merged_gfa,
            polished_record,
            &polished_fasta,
            &mapping,
            &node_classes,
            &blastn,
            &paths,
            options.keep_debug,
        )?;
        (
            sequence_plan,
            node_classes,
            vec![
                ("reconstruction_generator", "generator2_node_projection".to_string()),
                ("reconstruction_scenario", "no_recorded_resolve_path".to_string()),
                ("homologous_core_applied", "false".to_string()),
                ("homologous_core_trims", ".".to_string()),
                ("homologous_core_note", ".".to_string()),
            ],
        )
    };
    if let Some(correction) = &structural_corrections {
        for (key, value) in &mut path_summary {
            if *key == "topology_contract" {
                *value = "PASS_checked_lineage_with_explicit_structure_corrections".to_string();
            }
        }
        path_summary.extend([
            (
                "structure_corrections_applied",
                correction.count.to_string(),
            ),
            (
                "structure_removed_links",
                correction.removed_links.to_string(),
            ),
            ("structure_added_links", correction.added_links.to_string()),
            (
                "structure_reversed_node_sequences",
                correction.reversed_nodes.to_string(),
            ),
            (
                "structure_removed_nodes",
                correction.removed_nodes.to_string(),
            ),
            ("structure_original_walk", correction.old_walk.clone()),
            ("structure_corrected_walk", correction.new_walk.clone()),
            (
                "structure_evidence",
                "supplied_records_hash_verified;biological_support_reviewed_upstream".to_string(),
            ),
        ]);
    }
    let sequence_projection_method = sequence_plan.method.clone();
    if let Some(repeat_node) = options.repeat_resolution_node.as_deref() {
        if sequence_projection_method == "minimap2_fallback" {
            return Err(OrgraftError::InvalidArgument(format!(
                "repeat-constrained rebuild for `{repeat_node}` could not anchor every topology repeat to two polished FASTA copies; refusing the non-partitioning minimap2 fallback"
            )));
        }
    }
    let mapping = sequence_plan.mapping;
    let mut verified_gfa = make_verified_gfa_from_sequences(
        &merged_gfa,
        &sequence_plan.sequences,
        &mapping,
        &node_classes,
        &sequence_projection_method,
    );
    let validate_data_dir = discover_validate_data_dir(&polished_fasta);
    let coverage_rows =
        compute_node_remapped_coverage(&mapping, &node_classes, validate_data_dir.as_deref())?;
    annotate_gfa_with_remapped_coverage(&mut verified_gfa, &coverage_rows);
    // Spanning-read templates must come from assembly node sequences. Projecting
    // them through the currently tested linear FASTA would bake that FASTA's
    // repeat pairing into the evidence used to judge the pairing.
    let repeat_path_rows = repeat_path_support_rows(
        &support_source_gfa,
        &node_classes,
        reads.as_deref(),
        &minimap2,
        options.threads,
        &paths.debug_dir,
        options.keep_debug,
    )?;
    attach_repeat_path_support_paths(&mut verified_gfa, &repeat_path_rows);
    write_repeat_path_support(&paths.repeat_path_support, &repeat_path_rows)?;
    // Node coverage alone cannot prove that the oriented graph spells the
    // accepted sequence. In the no-contract path, check the final graph itself
    // before publishing either accepted output. Alternate repeat L edges need
    // not all occur in the one accepted circular walk.
    let mut graph_replay_status = None;
    if recorded_contract.is_none() {
        let candidate = verified_gfa.to_bytes()?;
        let replay = if verified_gfa.links.is_empty() {
            overlap_only::validate_linear_singleton(&verified_gfa, &polished_record.sequence)
        } else {
            overlap_only::validate_existing_circle(&candidate, &polished_record.sequence)
        };
        let replay_report = paths.logs_dir.join(format!(
            "rebuild_{}_graph_replay.tsv",
            options.subgraph
        ));
        match replay {
            Ok(proof) => {
                graph_replay_status = Some(proof.status);
                fs::write(
                    &replay_report,
                    format!(
                        "section\tkey\tvalue\nrun\tstatus\t{}\nproof\tpath_type\t{}\nproof\tused_segments\t{}\nproof\tused_links\t{}\nproof\tall_links_required\tfalse\nproof\tP_graph_adjacency\tPASS\nproof\tread_support\tnot_recomputed\nproof\tsame_method_resolve\t{}\nproof\tsearch_states\t{}\nproof\tsearch_starts\t{}\nproof\tsearch_transitions\t{}\n",
                        proof.status,
                        proof.path_type,
                        proof.used_segments,
                        proof.used_links,
                        if proof.status == "PASS" { "not_run" } else { "required" },
                        proof.states,
                        proof.starts,
                        proof.transitions,
                    ),
                )?;
            }
            Err(error) => {
                fs::write(
                    &replay_report,
                    format!(
                        "section\tkey\tvalue\nrun\tstatus\tREVIEW\nrun\treason\t{}\nproof\tread_support\tnot_recomputed\n",
                        error.to_string().replace(['\t', '\n'], " ")
                    ),
                )?;
                return Err(OrgraftError::InvalidArgument(format!(
                    "[REBUILD_SCENARIO=graph_replay_review] final graph cannot be accepted: {error}; no rebuilt GFA or FASTA written"
                )));
            }
        }
    }
    verified_gfa.write(&paths.verified_gfa)?;

    let repeat_resolution = if let Some(repeat_node) = options.repeat_resolution_node.as_deref() {
        let decision = select_dominant_repeat_pairing(&repeat_path_rows, repeat_node)?;
        let constrained = crate::commands::resolve::select_repeat_path_constrained_candidate(
            &crate::commands::resolve::RepeatPathConstraintRequest {
                component_gfa: &paths.verified_gfa,
                reference_fasta: &polished_fasta,
                repeat_node,
                pairing_signature: &decision.pairing_signature,
                max_states: 5_000,
                max_candidates: 100,
            },
        )?;
        let (sequence, orientation, rotation_offset, exact_position_matches) =
            normalize_repeat_constrained_sequence(
                &constrained.sequence,
                &polished_record.sequence,
                &constrained.boundary_offsets,
            )?;
        write_single_fasta(&paths.verified_fasta, &polished_record.header, &sequence)?;
        let local_candidate = constrained.selected_candidate.clone();
        let selected_candidate = source_repeat_candidate(
            prior_repeat_resolution.as_ref(),
            &decision,
            &local_candidate,
        );
        let outcome = RepeatResolutionOutcome {
            decision,
            selected_candidate,
            local_candidate,
            candidate_count: constrained.candidate_count,
            filtered_candidate_count: constrained.filtered_candidate_count,
            selection_method: constrained.selection_method,
            score: constrained.score,
            continuous_fraction: constrained.continuous_fraction,
            ordered_nodes: constrained.ordered_nodes,
            selected_order: constrained.selected_order,
            oriented_nodes: constrained.oriented_nodes,
            orientation,
            rotation_offset,
            exact_position_matches,
            candidates: constrained.candidates,
        };
        write_repeat_resolution(&paths.repeat_resolution, Some(&outcome))?;
        Some(outcome)
    } else {
        write_single_fasta(
            &paths.verified_fasta,
            &polished_record.header,
            &polished_record.sequence,
        )?;
        write_repeat_resolution(&paths.repeat_resolution, None)?;
        None
    };
    write_node_fasta(&paths.verified_nodes_fasta, &verified_gfa)?;
    let (gfa_editor_cli, image_exports) =
        if let Some(image_reference_fasta) = image_reference_fasta.as_deref() {
            match resolve_gfa_editor_cli(&soft_paths) {
                Ok(gfa_editor_cli) => {
                    let image_exports = export_gfa_reference_images(
                        &paths,
                        image_reference_fasta,
                        &gfa_editor_cli,
                        &soft_paths,
                    );
                    for row in image_exports.iter().filter(|row| row.status != "written") {
                        eprintln!(
                            "Warning: optional GFA_Editor {} export {} for {}; see {}",
                            row.format,
                            row.status.replace('_', " "),
                            paths.verified_gfa.display(),
                            paths.run_report.display()
                        );
                    }
                    (Some(gfa_editor_cli), image_exports)
                }
                Err(error) => {
                    eprintln!("Warning: optional GFA_Editor image export skipped: {error}");
                    (
                        None,
                        skipped_gfa_reference_images(&paths, &error.to_string()),
                    )
                }
            }
        } else {
            (None, Vec::new())
        };
    let consistency = CoordinateConsistency::new(
        polished_record,
        &raw_gfa,
        &merged_gfa,
        &verified_gfa,
        &mapping,
    );

    let raw_stats = graph_stats(&raw_gfa);
    let merged_stats = graph_stats(&merged_gfa);
    let verified_stats = graph_stats(&verified_gfa);
    let status = if graph_replay_status == Some("REVIEW_linear_no_closure") {
        "WARN"
    } else if recorded_contract.is_some()
        || (consistency.unmapped_segments == 0 && consistency.linear_tiling_status() == "PASS")
    {
        "PASS"
    } else {
        "WARN"
    };
    let mut summary = vec![
        ("status", status.to_string()),
        (
            "graph_replay_status",
            graph_replay_status.unwrap_or("recorded_contract").to_string(),
        ),
        ("organelle", options.organelle.clone()),
        ("subgraph", options.subgraph.clone()),
        ("edited_gfa", edited_gfa.display().to_string()),
        ("polished_fasta", polished_fasta.display().to_string()),
        ("merge_mode", merge_mode.clone()),
        (
            "repeat_resolution",
            if repeat_resolution.is_some() {
                "spanning_read_constrained"
            } else {
                "disabled"
            }
            .to_string(),
        ),
        (
            "repeat_resolution_node",
            repeat_resolution
                .as_ref()
                .map(|outcome| outcome.decision.repeat_node.clone())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "repeat_resolution_candidate",
            repeat_resolution
                .as_ref()
                .map(|outcome| outcome.selected_candidate.clone())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "repeat_resolution_local_candidate",
            repeat_resolution
                .as_ref()
                .map(|outcome| outcome.local_candidate.clone())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "sequence_projection_method",
            sequence_projection_method.clone(),
        ),
        ("raw_segments", raw_stats.segments.to_string()),
        ("raw_links", raw_stats.links.to_string()),
        ("raw_bp", raw_stats.bp.to_string()),
        ("merged_segments", merged_stats.segments.to_string()),
        ("merged_links", merged_stats.links.to_string()),
        ("merged_bp", merged_stats.bp.to_string()),
        ("verified_segments", verified_stats.segments.to_string()),
        ("verified_links", verified_stats.links.to_string()),
        ("verified_bp", verified_stats.bp.to_string()),
        ("mapped_segments", consistency.mapped_segments.to_string()),
        (
            "unmapped_segments",
            consistency.unmapped_segments.to_string(),
        ),
        ("polished_length", consistency.polished_length.to_string()),
        ("covered_bases", consistency.covered_bases.to_string()),
        ("gap_bases", consistency.gap_bases.to_string()),
        (
            "multi_covered_bases",
            consistency.multi_covered_bases.to_string(),
        ),
        (
            "linear_tiling",
            consistency.linear_tiling_status().to_string(),
        ),
        (
            "coverage_fraction",
            format!("{:.6}", consistency.coverage_fraction),
        ),
        (
            "runtime_seconds",
            format!("{:.6}", started.elapsed().as_secs_f64()),
        ),
    ];
    summary.extend(path_summary);
    write_rebuild_extract(
        &paths.extract_report,
        &options.subgraph,
        &verified_gfa,
        &mapping,
        &node_classes,
        &sequence_plan.node_source_rows,
        &coverage_rows,
    )?;
    write_run_report(
        &paths.run_report,
        &options.out_dir,
        &paths,
        options,
        &edited_gfa,
        &polished_fasta,
        reads.as_deref(),
        image_reference_fasta.as_deref(),
        validate_data_dir.as_deref(),
        &merge_mode,
        repeat_resolution.as_ref(),
        &sequence_projection_method,
        &minimap2,
        &blastn,
        gfa_editor_cli.as_deref(),
        &mapping,
        &image_exports,
        &started,
    )?;
    write_result_stats(
        &paths.result_stats,
        &options.subgraph,
        &summary,
        &raw_stats,
        &merged_stats,
        &verified_stats,
        &consistency,
        &coverage_rows,
        &repeat_path_rows,
    )?;
    if let Some(correction) = &structural_corrections {
        fs::write(
            paths.logs_dir.join(format!(
                "rebuild_{}_structure_corrections.tsv",
                options.subgraph
            )),
            &correction.record,
        )?;
    }

    println!("Wrote {}", paths.verified_gfa.display());
    println!("Wrote {}", paths.run_report.display());
    Ok(())
}

fn contract() -> CommandContract {
    CommandContract {
        command: "rebuild",
        origin: "high-quality graph generation after orgraft resolve and polish",
        purpose: "rebuild final verified graph and compact reports",
        inputs: &[
            "--organelle NAME",
            "--subgraph ID",
            "--edited-gfa FILE",
            "--polished-fasta FILE",
            "optional paired --checked-gfa FILE and --resolve-details FILE for checked topology and recorded-path sequence repair",
            "optional --ordinary-source-dir DIR replays a validated ordinary source transition only after recorded-path reconstruction fails",
            "optional paired --joint-ir-trial-gfa FILE and --joint-ir-evidence FILE for a fresh, hash-bound joint IR candidate using the original checked path",
            "optional --structure-corrections FILE binds graph/path/final FASTA hashes and explicit source evidence for reviewed structure corrections",
            "optional --reads FILE for complete flank-repeat-flank path support",
            "optional --repeat-resolution-node NODE for checkpoint-triggered path constraint",
            "optional --image-reference-fasta FILE",
            "optional --merged-gfa-template FILE",
            "soft_paths.txt containing minimap2, blastn, and optional gfa_editor_cli",
        ],
        outputs: &[
            "subgraph_001/rebuild_subgraph_001.gfa",
            "subgraph_001/rebuild_subgraph_001.fasta",
            "subgraph_001/rebuild_subgraph_001_nodes.fasta",
            "subgraph_001/rebuild_subgraph_001.pdf",
            "subgraph_001/rebuild_subgraph_001.svg",
            "logs/rebuild_subgraph_001_extract.tsv",
            "logs/rebuild_subgraph_001_run_report.tsv",
            "logs/rebuild_subgraph_001_result_stats.tsv",
            "logs/rebuild_subgraph_001_repeat_path_support.tsv",
            "logs/rebuild_subgraph_001_repeat_resolution.tsv",
        ],
        notes: &[
            "rebuild logic is Rust; optional PDF/SVG export calls the GFA_Editor CLI",
            "paired checked graph and resolve details keep the recorded path and shared nodes fixed, restoring only checked-supported topology; unresolved constraints stop the rebuild",
            "recorded-path mode proves the graph spell against the final FASTA; coordinate tiling WARN remains diagnostic and exact closing overlaps are represented once",
            "auto merge follows the conservative non-repeat linear compaction used by resolve/verified_gfa",
            "repeat resolution is disabled unless one checkpoint-triggered node is constrained by complete spanning-read support",
            "minimap2 projects merged graph nodes; blastn self-BLAST finds maximal exact repeat cores",
        ],
    }
}

#[derive(Debug, Clone)]
struct RebuildOptions {
    edited_gfa: PathBuf,
    polished_fasta: PathBuf,
    checked_gfa: Option<PathBuf>,
    resolve_details: Option<PathBuf>,
    ordinary_source_dir: Option<PathBuf>,
    joint_ir_trial_gfa: Option<PathBuf>,
    joint_ir_evidence: Option<PathBuf>,
    structure_corrections: Option<PathBuf>,
    reads: Option<PathBuf>,
    out_dir: PathBuf,
    soft_paths: PathBuf,
    image_reference_fasta: Option<PathBuf>,
    merged_gfa_template: Option<PathBuf>,
    minimap2: Option<PathBuf>,
    blastn: Option<PathBuf>,
    repeat_resolution_node: Option<String>,
    organelle: String,
    subgraph: String,
    threads: usize,
    force: bool,
    keep_debug: bool,
    overlap_only: bool,
}

impl RebuildOptions {
    fn from_args(args: &[String]) -> Result<Self, OrgraftError> {
        let mut edited_gfa = None;
        let mut polished_fasta = None;
        let mut checked_gfa = None;
        let mut resolve_details = None;
        let mut ordinary_source_dir = None;
        let mut joint_ir_trial_gfa = None;
        let mut joint_ir_evidence = None;
        // Retained but always None while --structure-corrections is disabled (see parser).
        let structure_corrections: Option<PathBuf> = None;
        let mut reads = None;
        let mut out_dir = PathBuf::from(DEFAULT_OUT_DIR);
        let mut soft_paths = PathBuf::from(DEFAULT_SOFT_PATHS);
        let mut image_reference_fasta = None;
        let mut merged_gfa_template = None;
        let mut minimap2 = None;
        let mut blastn = None;
        let mut repeat_resolution_node = None;
        let mut organelle = DEFAULT_ORGANELLE.to_string();
        let mut subgraph = DEFAULT_SUBGRAPH.to_string();
        let mut threads = DEFAULT_THREADS;
        let mut force = false;
        let mut keep_debug = false;
        let mut overlap_only = false;
        let mut out_dir_given = false;

        let mut index = 0usize;
        while index < args.len() {
            let arg = &args[index];
            match arg.as_str() {
                "--edited-gfa" | "--checked-draft-gfa" | "--draft-gfa" | "--raw-gfa" | "--gfa" => {
                    edited_gfa = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--polished-fasta" | "--verified-fasta" | "--final-fasta" | "--fasta" => {
                    polished_fasta = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--reads" => {
                    reads = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--checked-gfa" => {
                    checked_gfa = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--ordinary-source-dir" => {
                    ordinary_source_dir = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--resolve-details" => {
                    resolve_details = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--joint-ir-trial-gfa" => {
                    joint_ir_trial_gfa = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--joint-ir-evidence" => {
                    joint_ir_evidence = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--structure-corrections" => {
                    // TEMPORARILY DISABLED (decision 2026-09-12): manual, evidence-backed
                    // graph editing (link/reverse/remove/walk/node_sequence/split_node) is
                    // retained in rebuild_structure.rs but is NOT open in normal runs — a
                    // rebuild must not hand-edit the source graph from other-software
                    // results. Unreplayable candidates are flagged for review without
                    // assigning a source-graph cause (see reconstruction_failure_tag).
                    // Re-enable by restoring value parsing here and the apply block in run().
                    let _path = required_value(args, &mut index, arg)?;
                    return Err(OrgraftError::InvalidArgument(
                        "--structure-corrections is temporarily disabled: source-graph structural errors are detection-only".to_string(),
                    ));
                }
                "--out-dir" | "--output-dir" => {
                    out_dir = PathBuf::from(required_value(args, &mut index, arg)?);
                    out_dir_given = true;
                }
                "--overlap-only" => overlap_only = true,
                "--soft-paths" => {
                    soft_paths = PathBuf::from(required_value(args, &mut index, arg)?);
                }
                "--image-reference-fasta" | "--reference-fasta" => {
                    image_reference_fasta =
                        Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--merged-gfa-template" | "--merged-template" => {
                    merged_gfa_template =
                        Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--blastn" => {
                    blastn = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--minimap2" => {
                    minimap2 = Some(PathBuf::from(required_value(args, &mut index, arg)?));
                }
                "--repeat-resolution-node" => {
                    repeat_resolution_node =
                        Some(required_value(args, &mut index, arg)?.to_string());
                }
                "--organelle" => {
                    organelle = required_value(args, &mut index, arg)?.to_string();
                }
                "--subgraph" | "--ring" => {
                    subgraph = parse_subgraph_id(required_value(args, &mut index, arg)?)?;
                }
                "--threads" => {
                    threads = parse_usize(required_value(args, &mut index, arg)?, arg)?;
                    if threads == 0 {
                        return Err(OrgraftError::InvalidArgument(
                            "--threads must be at least 1".to_string(),
                        ));
                    }
                }
                "--force" => force = true,
                "--keep-debug" => keep_debug = true,
                other => {
                    return Err(OrgraftError::InvalidArgument(format!(
                        "unknown orgraft rebuild option `{other}`"
                    )));
                }
            }
            index += 1;
        }

        if ordinary_source_dir.is_some() && (checked_gfa.is_none() || reads.is_none() || force || overlap_only || joint_ir_trial_gfa.is_some()) {
            return Err(OrgraftError::InvalidArgument("ordinary source transition requires checked/resolve/reads and fresh normal output".into()));
        }
        if structure_corrections.is_some() && checked_gfa.is_none() {
            return Err(OrgraftError::InvalidArgument(
                "--structure-corrections requires --checked-gfa and --resolve-details".to_string(),
            ));
        }
        if checked_gfa.is_some() != resolve_details.is_some() {
            return Err(OrgraftError::InvalidArgument(
                "--checked-gfa and --resolve-details must be supplied together".to_string(),
            ));
        }
        if joint_ir_trial_gfa.is_some() != joint_ir_evidence.is_some()
            || (joint_ir_trial_gfa.is_some()
                && (checked_gfa.is_none() || force || overlap_only))
        {
            return Err(OrgraftError::InvalidArgument(
                "joint IR rebuild requires both trial GFA and evidence with original checked graph and resolve details, and a fresh output".to_string(),
            ));
        }
        if checked_gfa.is_some()
            && (repeat_resolution_node.is_some() || merged_gfa_template.is_some())
        {
            return Err(OrgraftError::InvalidArgument(
                "recorded-path rebuild cannot be combined with --repeat-resolution-node or --merged-gfa-template".to_string(),
            ));
        }
        if overlap_only
            && (!out_dir_given
                || force
                || checked_gfa.is_some()
                || resolve_details.is_some()
                || structure_corrections.is_some()
                || merged_gfa_template.is_some()
                || repeat_resolution_node.is_some()
                || reads.is_some()
                || image_reference_fasta.is_some())
        {
            return Err(OrgraftError::InvalidArgument(
                "--overlap-only requires an explicit fresh --out-dir and final --edited-gfa/--polished-fasta; it cannot combine with --force, graph projection, reads, or other rebuild modes".to_string(),
            ));
        }
        Ok(Self {
            edited_gfa: edited_gfa.ok_or_else(|| {
                OrgraftError::InvalidArgument("missing --edited-gfa FILE".to_string())
            })?,
            polished_fasta: polished_fasta.ok_or_else(|| {
                OrgraftError::InvalidArgument("missing --polished-fasta FILE".to_string())
            })?,
            checked_gfa,
            resolve_details,
            ordinary_source_dir,
            joint_ir_trial_gfa,
            joint_ir_evidence,
            structure_corrections,
            reads,
            out_dir,
            soft_paths,
            image_reference_fasta,
            merged_gfa_template,
            minimap2,
            blastn,
            repeat_resolution_node,
            organelle,
            subgraph,
            threads,
            force,
            keep_debug,
            overlap_only,
        })
    }
}

#[derive(Debug, Clone)]
struct OutputPaths {
    subgraph_dir: PathBuf,
    logs_dir: PathBuf,
    debug_dir: PathBuf,
    verified_gfa: PathBuf,
    verified_fasta: PathBuf,
    verified_nodes_fasta: PathBuf,
    verified_pdf: PathBuf,
    verified_svg: PathBuf,
    extract_report: PathBuf,
    run_report: PathBuf,
    result_stats: PathBuf,
    repeat_path_support: PathBuf,
    repeat_resolution: PathBuf,
}

impl OutputPaths {
    fn new(root: &Path, subgraph: &str) -> Self {
        let stem = format!("rebuild_{subgraph}");
        let subgraph_dir = root.join(subgraph);
        let logs_dir = root.join("logs");
        Self {
            subgraph_dir: subgraph_dir.clone(),
            logs_dir: logs_dir.clone(),
            debug_dir: root.join("debug").join(subgraph),
            verified_gfa: subgraph_dir.join(format!("{stem}.gfa")),
            verified_fasta: subgraph_dir.join(format!("{stem}.fasta")),
            verified_nodes_fasta: subgraph_dir.join(format!("{stem}_nodes.fasta")),
            verified_pdf: subgraph_dir.join(format!("{stem}.pdf")),
            verified_svg: subgraph_dir.join(format!("{stem}.svg")),
            extract_report: logs_dir.join(format!("{stem}_extract.tsv")),
            run_report: logs_dir.join(format!("{stem}_run_report.tsv")),
            result_stats: logs_dir.join(format!("{stem}_result_stats.tsv")),
            repeat_path_support: logs_dir.join(format!("{stem}_repeat_path_support.tsv")),
            repeat_resolution: logs_dir.join(format!("{stem}_repeat_resolution.tsv")),
        }
    }

    fn owned_paths(&self, subgraph: &str) -> Vec<PathBuf> {
        let mut paths = vec![self.subgraph_dir.clone(), self.debug_dir.clone()];
        // Match complete filenames, not a prefix: component IDs may themselves
        // contain underscores or be a prefix of a sibling component ID.
        for suffix in ["extract", "run_report", "result_stats", "repeat_path_support",
            "repeat_resolution", "graph_replay", "structure_corrections", "overlap_only", "ordinary_source"] {
            paths.push(self.logs_dir.join(format!("rebuild_{subgraph}_{suffix}.tsv")));
        }
        paths
    }

    fn prepare(&self, options: &RebuildOptions) -> Result<Option<PriorRepeatResolution>, OrgraftError> {
        if matches!(options.subgraph.as_str(), "logs" | "debug") {
            return Err(OrgraftError::InvalidArgument("--subgraph conflicts with a shared rebuild directory".into()));
        }
        for shared in [&options.out_dir, &self.logs_dir, self.debug_dir.parent().unwrap()] {
            let path = checked_rebuild_output_path(shared)?;
            if path.exists() && !path.is_dir() {
                return Err(OrgraftError::InvalidArgument(format!("rebuild output directory is a file: {}", path.display())));
            }
        }
        let owned = self.owned_paths(&options.subgraph);
        let destinations = owned.iter().map(|p| checked_rebuild_output_path(p)).collect::<Result<Vec<_>, _>>()?;
        let mut inputs = vec![options.edited_gfa.as_path(), options.polished_fasta.as_path(), options.soft_paths.as_path()];
        inputs.extend([
            options.checked_gfa.as_deref(), options.resolve_details.as_deref(), options.ordinary_source_dir.as_deref(),
            options.joint_ir_trial_gfa.as_deref(), options.joint_ir_evidence.as_deref(),
            options.structure_corrections.as_deref(), options.reads.as_deref(),
            options.image_reference_fasta.as_deref(), options.merged_gfa_template.as_deref(),
            options.minimap2.as_deref(), options.blastn.as_deref(),
        ].into_iter().flatten());
        for input in inputs {
            let lexical = absolute_rebuild_path(input)?;
            let canonical = if input.exists() { Some(fs::canonicalize(input)?) } else { None };
            for source in std::iter::once(&lexical).chain(canonical.iter()) {
                for (index, destination) in destinations.iter().enumerate() {
                    let overlaps = if index < 2 {
                        source.starts_with(destination) || destination.starts_with(source)
                    } else {
                        source == destination || (source.is_dir() && destination.starts_with(source))
                    };
                    if overlaps {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "rebuild output {} overlaps source {}; choose a separate component output",
                            destination.display(), input.display()
                        )));
                    }
                }
            }
        }
        let collisions = owned.iter().filter(|path| path.exists()).collect::<Vec<_>>();
        if !options.force {
            if let Some(path) = collisions.first() {
                return Err(OrgraftError::InvalidArgument(format!(
                    "{} already exists; use --force to replace only this component", path.display()
                )));
            }
        }
        // Preserve constrained-rebuild provenance before replacing this component.
        let prior = options.repeat_resolution_node.as_deref()
            .and_then(|node| read_prior_repeat_resolution(self, node));
        if options.force {
            for path in collisions { remove_path_if_exists(path)?; }
        }
        self.create(options.keep_debug)?;
        Ok(prior)
    }

    fn create(&self, keep_debug: bool) -> Result<(), OrgraftError> {
        fs::create_dir_all(self.subgraph_dir.parent().unwrap())?;
        fs::create_dir_all(&self.logs_dir)?;
        // Reserve the component before any native tool writes; a second process
        // cannot silently share the same component work directory.
        fs::create_dir(&self.subgraph_dir)?;
        if keep_debug {
            fs::create_dir_all(&self.debug_dir)?;
        }
        Ok(())
    }
}

fn absolute_rebuild_path(path: &Path) -> Result<PathBuf, OrgraftError> {
    use std::path::Component;
    let absolute = if path.is_absolute() { path.to_path_buf() } else { std::env::current_dir()?.join(path) };
    let mut normalized = PathBuf::new();
    for part in absolute.components() {
        match part {
            Component::CurDir => {},
            Component::ParentDir => { normalized.pop(); },
            _ => normalized.push(part.as_os_str()),
        }
    }
    Ok(normalized)
}

fn checked_rebuild_output_path(path: &Path) -> Result<PathBuf, OrgraftError> {
    // Inspect the original spelling too, so a symlink followed by `..` cannot
    // hide an output alias during lexical normalization.
    let absolute = if path.is_absolute() { path.to_path_buf() } else { std::env::current_dir()?.join(path) };
    let mut prefix = PathBuf::new();
    for part in absolute.components() {
        prefix.push(part.as_os_str());
        match fs::symlink_metadata(&prefix) {
            Ok(meta) if meta.file_type().is_symlink() => return Err(OrgraftError::InvalidArgument(format!(
                "rebuild output must not use a symlink: {}", prefix.display()
            ))),
            Ok(_) => {},
            Err(error) if error.kind() == io::ErrorKind::NotFound => {},
            Err(error) => return Err(error.into()),
        }
    }
    absolute_rebuild_path(path)
}

fn relative_output_path(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .to_string_lossy()
        .replace(std::path::MAIN_SEPARATOR, "/")
}

#[derive(Debug, Clone)]
// R4 note: the following GFA/FASTA types are exposed pub(crate) so the sibling
// sample-level `refine` module can reuse the exact same parser/writer
// and graph primitives. Visibility-only; no behaviour change (rebuild frozen).
pub(crate) struct FastaRecord {
    pub(crate) header: String,
    pub(crate) id: String,
    pub(crate) sequence: String,
}

#[derive(Debug, Clone)]
pub(crate) struct Segment {
    pub(crate) name: String,
    pub(crate) sequence: String,
    pub(crate) tags: Vec<String>,
}

#[derive(Debug, Clone)]
pub(crate) struct Link {
    pub(crate) from_name: String,
    pub(crate) from_orient: char,
    pub(crate) to_name: String,
    pub(crate) to_orient: char,
    pub(crate) overlap: String,
    pub(crate) tags: Vec<String>,
}

#[derive(Debug, Clone)]
pub(crate) struct Gfa {
    pub(crate) headers: Vec<Vec<String>>,
    pub(crate) segments: BTreeMap<String, Segment>,
    pub(crate) order: Vec<String>,
    pub(crate) links: Vec<Link>,
    pub(crate) other_lines: Vec<Vec<String>>,
}

impl Gfa {
    pub(crate) fn read(path: &Path) -> Result<Self, OrgraftError> {
        let text = fs::read_to_string(path).map_err(|error| {
            OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
        })?;
        let mut headers = Vec::new();
        let mut segments = BTreeMap::new();
        let mut order = Vec::new();
        let mut links = Vec::new();
        let mut other_lines = Vec::new();
        for (line_index, line) in text.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<String> = line.split('\t').map(str::to_string).collect();
            match fields.first().map(String::as_str) {
                Some("H") => headers.push(fields),
                Some("S") => {
                    if fields.len() < 3 {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "{}:{} malformed S line",
                            path.display(),
                            line_index + 1
                        )));
                    }
                    let name = fields[1].clone();
                    order.push(name.clone());
                    segments.insert(
                        name.clone(),
                        Segment {
                            name,
                            sequence: fields[2].clone(),
                            tags: fields[3..].to_vec(),
                        },
                    );
                }
                Some("L") => {
                    if fields.len() < 6 {
                        return Err(OrgraftError::InvalidArgument(format!(
                            "{}:{} malformed L line",
                            path.display(),
                            line_index + 1
                        )));
                    }
                    links.push(Link {
                        from_name: fields[1].clone(),
                        from_orient: parse_orient(&fields[2], path, line_index + 1)?,
                        to_name: fields[3].clone(),
                        to_orient: parse_orient(&fields[4], path, line_index + 1)?,
                        overlap: fields[5].clone(),
                        tags: fields[6..].to_vec(),
                    });
                }
                Some(_) => other_lines.push(fields),
                None => {}
            }
        }
        if segments.is_empty() {
            return Err(OrgraftError::InvalidArgument(format!(
                "{} contains no segments",
                path.display()
            )));
        }
        Ok(Self {
            headers,
            segments,
            order,
            links,
            other_lines,
        })
    }

    fn write_into<W: Write>(&self, out: &mut W) -> Result<(), OrgraftError> {
        for fields in &self.headers {
            writeln!(out, "{}", fields.join("\t"))?;
        }
        for name in &self.order {
            let segment = &self.segments[name];
            let mut fields = vec![
                "S".to_string(),
                segment.name.clone(),
                segment.sequence.clone(),
            ];
            fields.extend(segment.tags.clone());
            writeln!(out, "{}", fields.join("\t"))?;
        }
        for link in &self.links {
            let mut fields = vec![
                "L".to_string(),
                link.from_name.clone(),
                link.from_orient.to_string(),
                link.to_name.clone(),
                link.to_orient.to_string(),
                link.overlap.clone(),
            ];
            fields.extend(link.tags.clone());
            writeln!(out, "{}", fields.join("\t"))?;
        }
        for fields in &self.other_lines {
            writeln!(out, "{}", fields.join("\t"))?;
        }
        Ok(())
    }

    fn to_bytes(&self) -> Result<Vec<u8>, OrgraftError> {
        let mut out = Vec::new();
        self.write_into(&mut out)?;
        Ok(out)
    }

    pub(crate) fn write(&self, path: &Path) -> Result<(), OrgraftError> {
        let mut out = File::create(path)?;
        self.write_into(&mut out)
    }

    pub(crate) fn degrees(&self) -> BTreeMap<String, usize> {
        let mut degrees = BTreeMap::new();
        for name in &self.order {
            degrees.insert(name.clone(), 0usize);
        }
        for link in &self.links {
            *degrees.entry(link.from_name.clone()).or_insert(0) += 1;
            *degrees.entry(link.to_name.clone()).or_insert(0) += 1;
        }
        degrees
    }
}

#[derive(Debug, Clone)]
struct PathInfo {
    path: Vec<String>,
    orientations: Vec<char>,
    index: usize,
}

pub(crate) fn merge_unambiguous_gfa(raw: &Gfa) -> (Gfa, String) {
    let protected = repeat_like_nodes(raw);
    let mut adjacency: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for name in &raw.order {
        adjacency.insert(name.clone(), Vec::new());
    }
    for link in &raw.links {
        if protected.contains(&link.from_name) || protected.contains(&link.to_name) {
            continue;
        }
        adjacency
            .entry(link.from_name.clone())
            .or_default()
            .push(link.to_name.clone());
        adjacency
            .entry(link.to_name.clone())
            .or_default()
            .push(link.from_name.clone());
    }

    let mut visited = HashSet::new();
    let mut components = Vec::new();
    for name in &raw.order {
        if visited.contains(name) || protected.contains(name) {
            continue;
        }
        let mut stack = vec![name.clone()];
        let mut component = Vec::new();
        visited.insert(name.clone());
        while let Some(current) = stack.pop() {
            component.push(current.clone());
            for next in adjacency.get(&current).cloned().unwrap_or_default() {
                if visited.insert(next.clone()) {
                    stack.push(next);
                }
            }
        }
        components.push(component);
    }

    let mut path_records: Vec<(Vec<String>, Vec<char>)> = Vec::new();
    let mut used = HashSet::new();
    for component in components {
        if let Some(order) = ordered_linear_component(&component, &adjacency) {
            let orientations = orient_path(raw, &order);
            used.extend(order.iter().cloned());
            path_records.push((order, orientations));
        } else {
            let mut names = component;
            names.sort_by(|a, b| natural_cmp(a, b));
            for name in names {
                used.insert(name.clone());
                path_records.push((vec![name], vec!['+']));
            }
        }
    }
    for name in &raw.order {
        if used.insert(name.clone()) {
            path_records.push((vec![name.clone()], vec!['+']));
        }
    }

    let mut old_to_new = HashMap::new();
    let mut old_to_path_info = HashMap::new();
    let mut segments = BTreeMap::new();
    let mut order = Vec::new();
    let mut mergeable_path_count = 0usize;
    for (path, orientations) in &path_records {
        if path.len() > 1 {
            mergeable_path_count += 1;
        }
        let name = path.join("_");
        let mut sequence = String::new();
        for (index, (old_name, orient)) in path.iter().zip(orientations.iter()).enumerate() {
            let mut part = raw.segments[old_name].sequence.clone();
            if *orient == '-' {
                part = reverse_complement(&part);
            }
            sequence.push_str(&part);
            old_to_new.insert(old_name.clone(), name.clone());
            old_to_path_info.insert(
                old_name.clone(),
                PathInfo {
                    path: path.clone(),
                    orientations: orientations.clone(),
                    index,
                },
            );
        }
        let tags = vec![
            format!("LN:i:{}", sequence.len()),
            if path.len() > 1 {
                "SC:Z:linear_compaction"
            } else {
                "SC:Z:preserved_node"
            }
            .to_string(),
            "RR:Z:disabled".to_string(),
        ];
        order.push(name.clone());
        segments.insert(
            name.clone(),
            Segment {
                name,
                sequence,
                tags,
            },
        );
    }

    let mut links = Vec::new();
    let mut seen = HashSet::new();
    for link in &raw.links {
        if let Some(converted) = convert_link(link, &old_to_new, &old_to_path_info) {
            let key = merged_link_key(&converted);
            if seen.insert(key) {
                links.push(converted);
            }
        }
    }

    let mode = if mergeable_path_count == 0 {
        "input_already_merged_or_no_linear_compaction"
    } else {
        "auto_non_repeat_linear_compaction"
    };
    (
        Gfa {
            headers: raw.headers.clone(),
            segments,
            order,
            links,
            other_lines: raw.other_lines.clone(),
        },
        mode.to_string(),
    )
}

fn merged_link_key(link: &Link) -> String {
    let forward = format!(
        "{}{}:{}{}:{}",
        link.from_name, link.from_orient, link.to_name, link.to_orient, link.overlap
    );
    if link.from_name != link.to_name {
        return forward;
    }
    let reverse = format!(
        "{}{}:{}{}:{}",
        link.to_name,
        flip_orient(link.to_orient),
        link.from_name,
        flip_orient(link.from_orient),
        link.overlap
    );
    let canonical = if forward <= reverse { forward } else { reverse };
    format!("{canonical}:{}", link.tags.join("\t"))
}

fn repeat_like_nodes(gfa: &Gfa) -> BTreeSet<String> {
    let degrees = gfa.degrees();
    let branch: Vec<usize> = gfa
        .order
        .iter()
        .filter_map(|name| degrees.get(name).copied())
        .filter(|degree| *degree >= 3)
        .collect();
    if branch.is_empty() {
        return BTreeSet::new();
    }
    let branch_floor = *branch.iter().min().unwrap();
    let has_higher = branch.iter().any(|degree| *degree > branch_floor);
    gfa.order
        .iter()
        .filter(|name| {
            let degree = degrees.get(*name).copied().unwrap_or(0);
            if has_higher {
                degree > branch_floor
            } else {
                degree >= 3
            }
        })
        .cloned()
        .collect()
}

fn ordered_linear_component(
    component: &[String],
    adjacency: &BTreeMap<String, Vec<String>>,
) -> Option<Vec<String>> {
    let component_set: HashSet<String> = component.iter().cloned().collect();
    let mut sub_degrees = BTreeMap::new();
    for name in component {
        let degree = adjacency
            .get(name)
            .unwrap_or(&Vec::new())
            .iter()
            .filter(|next| component_set.contains(*next))
            .count();
        if degree > 2 {
            return None;
        }
        sub_degrees.insert(name.clone(), degree);
    }
    let mut endpoints: Vec<String> = sub_degrees
        .iter()
        .filter(|(_, degree)| **degree <= 1)
        .map(|(name, _)| name.clone())
        .collect();
    if component.len() > 1 && endpoints.len() != 2 {
        return None;
    }
    endpoints.sort_by(|a, b| natural_cmp(a, b));
    let mut order = Vec::new();
    let mut previous: Option<String> = None;
    let mut current = endpoints
        .first()
        .cloned()
        .or_else(|| component.first().cloned())?;
    loop {
        order.push(current.clone());
        let mut candidates: Vec<String> = adjacency
            .get(&current)
            .unwrap_or(&Vec::new())
            .iter()
            .filter(|next| component_set.contains(*next))
            .filter(|next| previous.as_ref() != Some(*next))
            .cloned()
            .collect();
        candidates.sort_by(|a, b| natural_cmp(a, b));
        let Some(next) = candidates.first().cloned() else {
            break;
        };
        previous = Some(current);
        current = next;
    }
    if order.len() == component.len() {
        Some(order)
    } else {
        None
    }
}

fn orient_path(gfa: &Gfa, order: &[String]) -> Vec<char> {
    for start_orient in ['+', '-'] {
        let mut orientations = vec![start_orient];
        let mut ok = true;
        for pair in order.windows(2) {
            let current = &pair[0];
            let next = &pair[1];
            let Some(orient) = next_orient(gfa, current, *orientations.last().unwrap(), next)
            else {
                ok = false;
                break;
            };
            orientations.push(orient);
        }
        if ok {
            return orientations;
        }
    }
    vec!['+'; order.len()]
}

fn next_orient(gfa: &Gfa, current: &str, current_orient: char, next: &str) -> Option<char> {
    for link in &gfa.links {
        if link.from_name == current && link.to_name == next && link.from_orient == current_orient {
            return Some(link.to_orient);
        }
        if link.to_name == current
            && link.from_name == next
            && flip_orient(link.to_orient) == current_orient
        {
            return Some(flip_orient(link.from_orient));
        }
    }
    None
}

fn convert_link(
    link: &Link,
    old_to_new: &HashMap<String, String>,
    old_to_path_info: &HashMap<String, PathInfo>,
) -> Option<Link> {
    let merged_from = old_to_new.get(&link.from_name)?;
    let merged_to = old_to_new.get(&link.to_name)?;
    if merged_from == merged_to && link.from_name != link.to_name {
        return None;
    }
    let from_endpoint = old_from_endpoint(link.from_orient);
    let to_endpoint = old_to_endpoint(link.to_orient);
    let merged_from_endpoint =
        map_old_endpoint_to_merged(&link.from_name, from_endpoint, old_to_path_info)?;
    let merged_to_endpoint =
        map_old_endpoint_to_merged(&link.to_name, to_endpoint, old_to_path_info)?;
    Some(Link {
        from_name: merged_from.clone(),
        from_orient: orient_for_from_endpoint(merged_from_endpoint),
        to_name: merged_to.clone(),
        to_orient: orient_for_to_endpoint(merged_to_endpoint),
        overlap: link.overlap.clone(),
        tags: link.tags.clone(),
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SegmentEndpoint {
    Left,
    Right,
}

fn old_from_endpoint(orient: char) -> SegmentEndpoint {
    match orient {
        '+' => SegmentEndpoint::Right,
        '-' => SegmentEndpoint::Left,
        _ => SegmentEndpoint::Right,
    }
}

fn old_to_endpoint(orient: char) -> SegmentEndpoint {
    match orient {
        '+' => SegmentEndpoint::Left,
        '-' => SegmentEndpoint::Right,
        _ => SegmentEndpoint::Left,
    }
}

fn orient_for_from_endpoint(endpoint: SegmentEndpoint) -> char {
    match endpoint {
        SegmentEndpoint::Right => '+',
        SegmentEndpoint::Left => '-',
    }
}

fn orient_for_to_endpoint(endpoint: SegmentEndpoint) -> char {
    match endpoint {
        SegmentEndpoint::Left => '+',
        SegmentEndpoint::Right => '-',
    }
}

fn map_old_endpoint_to_merged(
    old_name: &str,
    endpoint: SegmentEndpoint,
    old_to_path_info: &HashMap<String, PathInfo>,
) -> Option<SegmentEndpoint> {
    let info = old_to_path_info.get(old_name)?;
    let path_orient = info.orientations[info.index];
    if info.path.len() == 1 {
        return Some(match (path_orient, endpoint) {
            ('+', SegmentEndpoint::Left) | ('-', SegmentEndpoint::Right) => SegmentEndpoint::Left,
            ('+', SegmentEndpoint::Right) | ('-', SegmentEndpoint::Left) => SegmentEndpoint::Right,
            _ => endpoint,
        });
    }
    let first_endpoint = match path_orient {
        '+' => SegmentEndpoint::Left,
        '-' => SegmentEndpoint::Right,
        _ => SegmentEndpoint::Left,
    };
    if info.index == 0 && endpoint == first_endpoint {
        return Some(SegmentEndpoint::Left);
    }
    let last_endpoint = match path_orient {
        '+' => SegmentEndpoint::Right,
        '-' => SegmentEndpoint::Left,
        _ => SegmentEndpoint::Right,
    };
    if info.index == info.path.len() - 1 && endpoint == last_endpoint {
        return Some(SegmentEndpoint::Right);
    }
    None
}

#[derive(Debug, Clone)]
struct PafHit {
    query_name: String,
    query_length: usize,
    query_start: usize,
    query_end: usize,
    strand: char,
    target_start: usize,
    target_end: usize,
    matches: usize,
    block_length: usize,
    mapq: usize,
    substitution_matches: usize,
    substitution_mismatches: usize,
    primary: bool,
    exact_flanks: bool,
}

impl PafHit {
    fn identity(&self) -> f64 {
        if self.block_length == 0 {
            0.0
        } else {
            self.matches as f64 / self.block_length as f64
        }
    }
    fn aligned_fraction(&self) -> f64 {
        if self.query_length == 0 {
            0.0
        } else {
            (self.query_end.saturating_sub(self.query_start)) as f64 / self.query_length as f64
        }
    }
    fn substitution_identity(&self) -> f64 {
        let total = self
            .substitution_matches
            .saturating_add(self.substitution_mismatches);
        if total == 0 {
            self.identity()
        } else {
            self.substitution_matches as f64 / total as f64
        }
    }
    fn folded(&self, reference_len: usize) -> FoldedInterval {
        FoldedInterval::new(self.target_start, self.target_end, reference_len)
    }
}

#[derive(Debug, Clone)]
struct FoldedInterval {
    // Keep every traversal of the circular reference. Endpoints alone cannot
    // distinguish a short overlap from a complete turn followed by that overlap.
    parts: Vec<(usize, usize)>,
}

impl FoldedInterval {
    fn new(start_zero: usize, end_zero: usize, reference_len: usize) -> Self {
        let mut remaining = end_zero.saturating_sub(start_zero);
        let mut parts = Vec::new();
        if reference_len == 0 {
            return Self { parts };
        }
        let mut position = start_zero % reference_len;
        while remaining > 0 {
            let length = remaining.min(reference_len - position);
            parts.push((position + 1, position + length));
            remaining -= length;
            position = 0;
        }
        Self { parts }
    }
    fn parts_text(&self, _reference_len: usize) -> String {
        intervals_text(&self.parts)
    }
}

#[derive(Debug, Clone)]
struct NodeMapping {
    selected: Option<PafHit>,
    accepted_hits: Vec<PafHit>,
    all_hits: Vec<PafHit>,
}

#[derive(Debug, Clone)]
struct Mapping {
    by_node: BTreeMap<String, NodeMapping>,
    reference_len: usize,
    command: String,
    stderr: String,
}

/// Detection-only scenario tag emitted when both generators fail. A failed replay
/// after node projection proves that this candidate graph cannot be accepted; it
/// does not identify whether the source topology, projection, or path choice is
/// responsible. Keep that distinction for review instead of attributing a source
/// graph error without independent evidence. Historical logs retain their old tag.
fn reconstruction_failure_tag(projection_error: &OrgraftError) -> &'static str {
    let message = projection_error.to_string();
    if message.contains("replay length differs")
        || message.contains("does not equal the authoritative FASTA")
    {
        "[REBUILD_SCENARIO=graph_replay_review] "
    } else {
        "[REBUILD_SCENARIO=unresolved_after_both_generators] "
    }
}

#[cfg(test)]
mod tests_reconstruction_failure_tag {
    use super::*;

    #[test]
    fn failed_replay_does_not_assert_a_source_graph_error() {
        for reason in [
            "candidate graph replay length differs from the authoritative FASTA",
            "candidate graph replay does not equal the authoritative FASTA in the allowed direction and origin",
        ] {
            let error = OrgraftError::InvalidArgument(reason.to_string());
            assert_eq!(reconstruction_failure_tag(&error), "[REBUILD_SCENARIO=graph_replay_review] ");
        }
    }

    #[test]
    fn other_generator_failures_keep_the_unresolved_label() {
        let error = OrgraftError::InvalidArgument("search limit reached".to_string());
        assert_eq!(reconstruction_failure_tag(&error), "[REBUILD_SCENARIO=unresolved_after_both_generators] ");
    }
}

const UNIQUE_PROJECTION_FLANK_BP: usize = 100;

fn usable_node_projection_hit(hit: &PafHit, is_repeat: bool) -> bool {
    if hit.aligned_fraction() < 0.05 {
        return false;
    }
    hit.identity() >= 0.97
        || (!is_repeat
            && hit.mapq >= 20
            && hit.aligned_fraction() >= 0.75
            && hit.substitution_identity() >= 0.99)
}

fn unique_circular_gap_projection_hit(hits: &[PafHit], reference_len: usize) -> Option<PafHit> {
    if hits.len() != 2 || reference_len == 0 {
        return None;
    }
    let mut pair = hits.iter().collect::<Vec<_>>();
    pair.sort_by_key(|hit| hit.target_start);
    let (first, second) = (pair[0], pair[1]);
    let full_and_strong = |hit: &PafHit| {
        hit.query_start == 0
            && hit.query_end == hit.query_length
            && hit.query_length >= UNIQUE_PROJECTION_FLANK_BP * 2
            && hit.mapq == 0
            && hit.substitution_identity() >= 0.99
            && hit.exact_flanks
    };
    if !full_and_strong(first)
        || !full_and_strong(second)
        || first.primary == second.primary
        || first.strand != second.strand
        || first.query_length != second.query_length
        || first.target_end > reference_len
        || second.target_end > reference_len.saturating_mul(2)
        || second.target_start != first.target_start.checked_add(reference_len)?
        || second.target_end != first.target_end.checked_add(reference_len)?
        || first.folded(reference_len).parts != second.folded(reference_len).parts
    {
        return None;
    }
    Some(first.clone())
}

fn run_node_projection(
    gfa: &Gfa,
    record: &FastaRecord,
    minimap2: &Path,
    threads: usize,
    paths: &OutputPaths,
    keep_debug: bool,
) -> Result<Mapping, OrgraftError> {
    let temp_dir = if keep_debug {
        paths.debug_dir.clone()
    } else {
        temp_work_dir("orgraft-rebuild")
    };
    fs::create_dir_all(&temp_dir)?;
    let query_fasta = temp_dir.join("merged_nodes.fasta");
    let target_fasta = temp_dir.join("polished_doubled.fasta");
    write_node_fasta_from_order(&query_fasta, gfa)?;
    write_single_fasta(
        &target_fasta,
        &format!("{}_doubled", record.id),
        &format!("{}{}", record.sequence, record.sequence),
    )?;

    let command_text = format!(
        "{} -x asm5 -c --eqx -N 20 -t {} {} {}",
        minimap2.display(),
        threads,
        target_fasta.display(),
        query_fasta.display()
    );
    let output = Command::new(minimap2)
        .args([
            "-x",
            "asm5",
            "-c",
            "--eqx",
            "-N",
            "20",
            "-t",
            &threads.to_string(),
        ])
        .arg(&target_fasta)
        .arg(&query_fasta)
        .output()
        .map_err(|error| {
            OrgraftError::InvalidArgument(format!("failed to run minimap2: {error}"))
        })?;
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    if keep_debug {
        fs::write(paths.debug_dir.join("minimap2.paf"), &stdout)?;
        fs::write(paths.debug_dir.join("minimap2.stderr.log"), &stderr)?;
    }
    if !output.status.success() {
        if !keep_debug {
            let _ = fs::remove_dir_all(&temp_dir);
        }
        return Err(OrgraftError::InvalidArgument(format!(
            "minimap2 failed: {stderr}"
        )));
    }

    let mut raw_hits: BTreeMap<String, Vec<PafHit>> = BTreeMap::new();
    for line in stdout.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let hit = parse_paf(line)?;
        raw_hits
            .entry(hit.query_name.clone())
            .or_default()
            .push(hit);
    }
    let mut by_node = BTreeMap::new();
    let repeat_nodes = repeat_like_nodes(gfa);
    for name in &gfa.order {
        let mut hits = raw_hits.remove(name).unwrap_or_default();
        let is_repeat = repeat_nodes.contains(name);
        let circular_gap_hit = if is_repeat || hits.iter().any(|hit| usable_node_projection_hit(hit, false)) {
            None
        } else {
            unique_circular_gap_projection_hit(&hits, record.sequence.len())
        };
        hits.retain(|hit| usable_node_projection_hit(hit, is_repeat));
        if let Some(hit) = &circular_gap_hit {
            hits.push(hit.clone());
        }
        sort_hits(&mut hits);
        let accepted_hits = if let Some(hit) = circular_gap_hit {
            vec![hit]
        } else {
            dedupe_hits(
                hits.iter()
                    .filter(|hit| hit.identity() >= 0.99 && hit.aligned_fraction() >= 0.75)
                    .cloned()
                    .collect(),
                record.sequence.len(),
            )
        };
        by_node.insert(
            name.clone(),
            NodeMapping {
                selected: hits.first().cloned(),
                accepted_hits,
                all_hits: hits,
            },
        );
    }
    if !keep_debug {
        let _ = fs::remove_dir_all(&temp_dir);
    }
    Ok(Mapping {
        by_node,
        reference_len: record.sequence.len(),
        command: command_text,
        stderr,
    })
}

fn classify_nodes(gfa: &Gfa, mapping: &Mapping) -> BTreeMap<String, String> {
    let repeats = repeat_like_nodes(gfa);
    let mut classes = BTreeMap::new();
    for name in &gfa.order {
        let mapped_repeat = mapping
            .by_node
            .get(name)
            .map(|node| node.accepted_hits.len() > 1)
            .unwrap_or(false);
        let class = if repeats.contains(name) || mapped_repeat {
            "repeat_node"
        } else {
            "single_copy_node"
        };
        classes.insert(name.clone(), class.to_string());
    }
    classes
}

#[derive(Debug, Clone)]
struct VerifiedNodeSourceRow {
    node: String,
    node_class: String,
    source: String,
    length: usize,
    parts: String,
    notes: String,
}

#[derive(Debug, Clone)]
struct SequencePlan {
    sequences: BTreeMap<String, String>,
    mapping: Mapping,
    node_source_rows: Vec<VerifiedNodeSourceRow>,
    method: String,
}

fn make_verified_gfa_from_sequences(
    merged: &Gfa,
    sequences: &BTreeMap<String, String>,
    mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    method: &str,
) -> Gfa {
    let mut segments = BTreeMap::new();
    for name in &merged.order {
        let segment = &merged.segments[name];
        let node_mapping = &mapping.by_node[name];
        let sequence = sequences
            .get(name)
            .cloned()
            .unwrap_or_else(|| segment.sequence.clone());
        let interval_text = if node_mapping.accepted_hits.is_empty() {
            ".".to_string()
        } else {
            node_mapping
                .accepted_hits
                .iter()
                .map(|hit| {
                    let interval = hit.folded(mapping.reference_len);
                    format!(
                        "{}:{}",
                        interval.parts_text(mapping.reference_len),
                        hit.strand
                    )
                })
                .collect::<Vec<_>>()
                .join(",")
        };
        let mut tags = replace_tags(
            &segment.tags,
            &[
                format!("VC:Z:{method}"),
                "VZ:Z:verified_fasta".to_string(),
                format!("VS:Z:{interval_text}"),
                format!("NC:Z:{}", node_classes[name]),
                format!("LN:i:{}", sequence.len()),
                format!("OL:i:{}", segment.sequence.len()),
                format!("VL:i:{}", sequence.len()),
            ],
        );
        if node_classes[name] == "repeat_node" {
            tags = replace_tag(&tags, "RT:Z:unresolved".to_string());
            tags = replace_tag(&tags, "PM:Z:repeat_path_support_pending".to_string());
        }
        segments.insert(
            name.clone(),
            Segment {
                name: name.clone(),
                sequence,
                tags,
            },
        );
    }
    Gfa {
        headers: merged.headers.clone(),
        segments,
        order: merged.order.clone(),
        links: merged.links.clone(),
        other_lines: merged.other_lines.clone(),
    }
}

#[derive(Debug, Clone)]
struct BlastRow {
    pident: f64,
    length: usize,
    mismatch: usize,
    gapopen: usize,
    qstart: usize,
    qend: usize,
    sstart: usize,
    send: usize,
    bitscore: f64,
}

#[derive(Debug, Clone)]
struct RepeatCopy {
    start: usize,
    end: usize,
    strand: char,
}

#[derive(Debug, Clone)]
struct RepeatCorePair {
    length: usize,
    source_blast_length: usize,
    source_pident: f64,
    source_gapopen: usize,
    source_mismatch: usize,
    boundary_extension_left: usize,
    boundary_extension_right: usize,
    copies: Vec<RepeatCopy>,
}

#[derive(Debug, Clone)]
struct RepeatGap {
    left_repeat_node: String,
    left_copy_index: usize,
    right_repeat_node: String,
    right_copy_index: usize,
    parts: Vec<(usize, usize)>,
}

fn infer_verified_sequences_from_repeat_cores(
    merged: &Gfa,
    record: &FastaRecord,
    verified_fasta: &Path,
    mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    blastn: &Path,
    paths: &OutputPaths,
    keep_debug: bool,
) -> Result<SequencePlan, OrgraftError> {
    let repeat_count = node_classes
        .values()
        .filter(|class| class.as_str() == "repeat_node")
        .count();
    if repeat_count == 0 {
        return Ok(fallback_sequence_plan(
            merged,
            record,
            mapping,
            node_classes,
            "minimap2_fallback",
        ));
    }

    let repeat_core_min_alignment_bp = repeat_core_min_alignment_bp(merged, node_classes);
    let core_pairs = infer_repeat_core_pairs(
        verified_fasta,
        record,
        blastn,
        repeat_count,
        repeat_core_min_alignment_bp,
    )?;
    if keep_debug {
        write_self_blast_repeat_core_debug(
            &paths.debug_dir.join("self_blast_repeat_cores.tsv"),
            &core_pairs,
        )?;
    }
    if core_pairs.len() < repeat_count {
        return Ok(fallback_sequence_plan(
            merged,
            record,
            mapping,
            node_classes,
            "minimap2_fallback",
        ));
    }
    let repeat_assignments =
        assign_repeat_nodes_to_core_pairs(merged, node_classes, mapping, &core_pairs);
    if repeat_assignments.len() < repeat_count {
        return Ok(fallback_sequence_plan(
            merged,
            record,
            mapping,
            node_classes,
            "minimap2_fallback",
        ));
    }

    let gaps = build_repeat_core_gaps(&repeat_assignments, record.sequence.len());
    let single_assignments = assign_single_copy_nodes_to_gaps(merged, node_classes, mapping, &gaps);
    let mut sequences = BTreeMap::new();
    let mut by_node = BTreeMap::new();
    let mut rows = Vec::new();

    for name in &merged.order {
        if node_classes[name] == "repeat_node" {
            let Some(pair) = repeat_assignments.get(name) else {
                continue;
            };
            let (copy, strand) = choose_repeat_copy_for_node(name, pair, mapping);
            let sequence = record_interval_sequence(record, copy.start, copy.end, strand);
            let synthetic_hits =
                synthetic_hits_for_repeat_copies(name, &sequence, record, &pair.copies);
            sequences.insert(name.clone(), sequence.clone());
            by_node.insert(
                name.clone(),
                NodeMapping {
                    selected: synthetic_hits.first().cloned(),
                    accepted_hits: synthetic_hits.clone(),
                    all_hits: synthetic_hits,
                },
            );
            rows.push(VerifiedNodeSourceRow {
                node: name.clone(),
                node_class: "repeat_node".to_string(),
                source: "verified_fasta_self_blast_repeat_core".to_string(),
                length: sequence.len(),
                parts: format!("{}-{}:{}", copy.start, copy.end, strand),
                notes: format!(
                    "core_copies={};blast_length={};blast_pident={:.3};boundary_extension_left={};boundary_extension_right={};repeat_pair_anchor={}",
                    pair.copies
                        .iter()
                        .map(|item| format!("{}-{}:{}", item.start, item.end, item.strand))
                        .collect::<Vec<_>>()
                        .join(","),
                    pair.source_blast_length,
                    pair.source_pident,
                    pair.boundary_extension_left,
                    pair.boundary_extension_right,
                    if mapping.by_node[name].accepted_hits.len() >= 2 {
                        "strict_block_identity"
                    } else {
                        "indel_tolerant_secondary"
                    }
                ),
            });
        }
    }

    for name in &merged.order {
        if node_classes[name] != "single_copy_node" {
            continue;
        }
        let Some(gap) = single_assignments.get(name) else {
            continue;
        };
        let strand = best_orientation_for_parts(name, &gap.parts, mapping);
        let sequence = verified_sequence_from_parts(record, &gap.parts, strand);
        let synthetic_hits =
            synthetic_hits_for_parts(name, sequence.len(), record, &gap.parts, strand);
        sequences.insert(name.clone(), sequence.clone());
        by_node.insert(
            name.clone(),
            NodeMapping {
                selected: synthetic_hits.first().cloned(),
                accepted_hits: synthetic_hits.clone(),
                all_hits: synthetic_hits,
            },
        );
        rows.push(VerifiedNodeSourceRow {
            node: name.clone(),
            node_class: "single_copy_node".to_string(),
            source: "verified_fasta_between_repeat_cores".to_string(),
            length: sequence.len(),
            parts: format!(
                "{}:{}",
                gap.parts
                    .iter()
                    .map(|(start, end)| format!("{start}-{end}"))
                    .collect::<Vec<_>>()
                    .join(","),
                strand
            ),
            notes: format!(
                "left={}:copy{};right={}:copy{}",
                gap.left_repeat_node,
                gap.left_copy_index,
                gap.right_repeat_node,
                gap.right_copy_index
            ),
        });
    }

    for name in &merged.order {
        if sequences.contains_key(name) {
            continue;
        }
        let segment = &merged.segments[name];
        let original = mapping.by_node.get(name);
        let hit = original.and_then(|node| node.selected.clone());
        let (sequence, accepted_hits, source, parts, notes) = if let Some(hit) = hit {
            let sequence = extract_verified_sequence(record, &hit);
            (
                sequence,
                vec![hit.clone()],
                "verified_fasta_minimap2_fallback".to_string(),
                format!("{}-{}:{}", hit.target_start + 1, hit.target_end, hit.strand),
                "fallback".to_string(),
            )
        } else {
            (
                segment.sequence.clone(),
                Vec::new(),
                "unmapped_original_sequence".to_string(),
                ".".to_string(),
                "fallback_no_alignment".to_string(),
            )
        };
        sequences.insert(name.clone(), sequence.clone());
        by_node.insert(
            name.clone(),
            NodeMapping {
                selected: accepted_hits.first().cloned(),
                accepted_hits: accepted_hits.clone(),
                all_hits: accepted_hits,
            },
        );
        rows.push(VerifiedNodeSourceRow {
            node: name.clone(),
            node_class: node_classes[name].clone(),
            source,
            length: sequence.len(),
            parts,
            notes,
        });
    }

    Ok(SequencePlan {
        sequences,
        mapping: Mapping {
            by_node,
            reference_len: record.sequence.len(),
            command: mapping.command.clone(),
            stderr: mapping.stderr.clone(),
        },
        node_source_rows: rows,
        method: "verified_fasta_self_blast_repeat_core".to_string(),
    })
}

fn fallback_sequence_plan(
    merged: &Gfa,
    record: &FastaRecord,
    mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    method: &str,
) -> SequencePlan {
    let mut sequences = BTreeMap::new();
    let mut rows = Vec::new();
    for name in &merged.order {
        let segment = &merged.segments[name];
        let node = &mapping.by_node[name];
        if let Some(hit) = &node.selected {
            let sequence = extract_verified_sequence(record, hit);
            sequences.insert(name.clone(), sequence.clone());
            rows.push(VerifiedNodeSourceRow {
                node: name.clone(),
                node_class: node_classes[name].clone(),
                source: "polished_fasta_minimap2".to_string(),
                length: sequence.len(),
                parts: format!("{}-{}:{}", hit.target_start + 1, hit.target_end, hit.strand),
                notes: format!(
                    "identity={:.6};aligned_fraction={:.6}",
                    hit.identity(),
                    hit.aligned_fraction()
                ),
            });
        } else {
            sequences.insert(name.clone(), segment.sequence.clone());
            rows.push(VerifiedNodeSourceRow {
                node: name.clone(),
                node_class: node_classes[name].clone(),
                source: "unmapped_original_sequence".to_string(),
                length: segment.sequence.len(),
                parts: ".".to_string(),
                notes: "selected_hit=missing".to_string(),
            });
        }
    }
    SequencePlan {
        sequences,
        mapping: mapping.clone(),
        node_source_rows: rows,
        method: method.to_string(),
    }
}

fn infer_repeat_core_pairs(
    verified_fasta: &Path,
    record: &FastaRecord,
    blastn: &Path,
    repeat_count: usize,
    min_alignment_bp: usize,
) -> Result<Vec<RepeatCorePair>, OrgraftError> {
    let mut rows = run_self_blastn_repeat_hits(verified_fasta, blastn)?;
    rows.retain(|row| row.length >= min_alignment_bp && row.pident >= 95.0);
    rows.sort_by(|a, b| {
        b.length
            .cmp(&a.length)
            .then_with(|| b.pident.partial_cmp(&a.pident).unwrap_or(Ordering::Equal))
            .then_with(|| {
                b.bitscore
                    .partial_cmp(&a.bitscore)
                    .unwrap_or(Ordering::Equal)
            })
    });

    let mut selected = Vec::new();
    let mut seen = HashSet::new();
    let mut used_intervals: Vec<(usize, usize)> = Vec::new();
    for row in rows {
        let key = core_pair_key(&row);
        if !seen.insert(key) {
            continue;
        }
        let (q_low, q_high) = sorted_pair(row.qstart, row.qend);
        let (s_low, s_high) = sorted_pair(row.sstart, row.send);
        if used_intervals.iter().any(|(old_start, old_end)| {
            interval_overlap(q_low, q_high, *old_start, *old_end) > 100
                || interval_overlap(s_low, s_high, *old_start, *old_end) > 100
        }) {
            continue;
        }
        let pair = repeat_core_pair_from_blast_row(&row, record);
        if pair.length == 0 {
            continue;
        }
        used_intervals.push((pair.copies[0].start, pair.copies[0].end));
        used_intervals.push((pair.copies[1].start, pair.copies[1].end));
        selected.push(pair);
        if selected.len() >= repeat_count {
            break;
        }
    }
    Ok(selected)
}

fn repeat_core_min_alignment_bp(merged: &Gfa, node_classes: &BTreeMap<String, String>) -> usize {
    let shortest_repeat_bp = merged
        .order
        .iter()
        .filter(|name| node_classes[*name] == "repeat_node")
        .filter_map(|name| merged.segments.get(name))
        .map(|segment| segment.sequence.len())
        .min()
        .unwrap_or(REPEAT_CORE_DEFAULT_ALIGNMENT_BP);
    repeat_core_alignment_threshold(shortest_repeat_bp)
}

fn repeat_core_alignment_threshold(shortest_repeat_bp: usize) -> usize {
    shortest_repeat_bp.saturating_mul(3).div_ceil(4).clamp(
        REPEAT_CORE_MIN_ALIGNMENT_BP,
        REPEAT_CORE_DEFAULT_ALIGNMENT_BP,
    )
}

fn run_self_blastn_repeat_hits(
    verified_fasta: &Path,
    blastn: &Path,
) -> Result<Vec<BlastRow>, OrgraftError> {
    let output = Command::new(blastn)
        .arg("-query")
        .arg(verified_fasta)
        .arg("-subject")
        .arg(verified_fasta)
        .arg("-outfmt")
        .arg("6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore")
        .arg("-dust")
        .arg("no")
        .arg("-soft_masking")
        .arg("false")
        .output()
        .map_err(|error| {
            OrgraftError::InvalidArgument(format!("failed to run blastn self alignment: {error}"))
        })?;
    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    if !output.status.success() {
        return Err(OrgraftError::InvalidArgument(format!(
            "blastn self alignment failed: {stderr}"
        )));
    }
    let mut rows = Vec::new();
    for line in stdout.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }
        let row = BlastRow {
            pident: fields[2].parse::<f64>().unwrap_or(0.0),
            length: fields[3].parse::<usize>().unwrap_or(0),
            mismatch: fields[4].parse::<usize>().unwrap_or(0),
            gapopen: fields[5].parse::<usize>().unwrap_or(0),
            qstart: fields[6].parse::<usize>().unwrap_or(0),
            qend: fields[7].parse::<usize>().unwrap_or(0),
            sstart: fields[8].parse::<usize>().unwrap_or(0),
            send: fields[9].parse::<usize>().unwrap_or(0),
            bitscore: fields[11].parse::<f64>().unwrap_or(0.0),
        };
        let (q_low, q_high) = sorted_pair(row.qstart, row.qend);
        let (s_low, s_high) = sorted_pair(row.sstart, row.send);
        if q_low == s_low && q_high == s_high {
            continue;
        }
        rows.push(row);
    }
    Ok(rows)
}

fn core_pair_key(row: &BlastRow) -> String {
    let (q_low, q_high) = sorted_pair(row.qstart, row.qend);
    let (s_low, s_high) = sorted_pair(row.sstart, row.send);
    let mut pairs = [(q_low, q_high), (s_low, s_high)];
    pairs.sort();
    format!(
        "{}-{}:{}-{}",
        pairs[0].0, pairs[0].1, pairs[1].0, pairs[1].1
    )
}

fn repeat_core_pair_from_blast_row(row: &BlastRow, record: &FastaRecord) -> RepeatCorePair {
    let (q_seq, q_strand, q_low, q_high) = blast_interval_sequence(record, row.qstart, row.qend);
    let (s_seq, s_strand, s_low, s_high) = blast_interval_sequence(record, row.sstart, row.send);
    let (q_offset, s_offset, core_length) = longest_exact_common_substring(&q_seq, &s_seq);
    let (q_start, q_end, q_core_strand) =
        map_oriented_offset_to_record(q_low, q_high, q_strand, q_offset, core_length);
    let (s_start, s_end, s_core_strand) =
        map_oriented_offset_to_record(s_low, s_high, s_strand, s_offset, core_length);
    let mut pair = RepeatCorePair {
        length: core_length,
        source_blast_length: row.length,
        source_pident: row.pident,
        source_gapopen: row.gapopen,
        source_mismatch: row.mismatch,
        boundary_extension_left: 0,
        boundary_extension_right: 0,
        copies: vec![
            RepeatCopy {
                start: q_start,
                end: q_end,
                strand: q_core_strand,
            },
            RepeatCopy {
                start: s_start,
                end: s_end,
                strand: s_core_strand,
            },
        ],
    };
    extend_repeat_core_pair_to_max_exact(&mut pair, record);
    pair
}

fn blast_interval_sequence(
    record: &FastaRecord,
    start: usize,
    end: usize,
) -> (String, char, usize, usize) {
    if start <= end {
        (
            record_interval_sequence(record, start, end, '+'),
            '+',
            start,
            end,
        )
    } else {
        (
            record_interval_sequence(record, end, start, '-'),
            '-',
            end,
            start,
        )
    }
}

fn longest_exact_common_substring(seq_a: &str, seq_b: &str) -> (usize, usize, usize) {
    let a = seq_a.as_bytes();
    let b = seq_b.as_bytes();
    let mut previous = vec![0usize; b.len() + 1];
    let mut best_start_a = 0usize;
    let mut best_start_b = 0usize;
    let mut best_len = 0usize;
    for (index_a, base_a) in a.iter().enumerate() {
        let mut current = vec![0usize; b.len() + 1];
        for (index_b, base_b) in b.iter().enumerate() {
            if base_a.eq_ignore_ascii_case(base_b) {
                current[index_b + 1] = previous[index_b] + 1;
                if current[index_b + 1] > best_len {
                    best_len = current[index_b + 1];
                    best_start_a = index_a + 1 - best_len;
                    best_start_b = index_b + 1 - best_len;
                }
            }
        }
        previous = current;
    }
    (best_start_a, best_start_b, best_len)
}

fn map_oriented_offset_to_record(
    start: usize,
    end: usize,
    strand: char,
    offset: usize,
    length: usize,
) -> (usize, usize, char) {
    if strand == '+' {
        let core_start = start + offset;
        (core_start, core_start + length.saturating_sub(1), '+')
    } else {
        let core_end = end.saturating_sub(offset);
        (core_end + 1 - length, core_end, '-')
    }
}

fn extend_repeat_core_pair_to_max_exact(pair: &mut RepeatCorePair, record: &FastaRecord) {
    while try_extend_repeat_core_pair(pair, record, "left") {
        pair.boundary_extension_left += 1;
    }
    while try_extend_repeat_core_pair(pair, record, "right") {
        pair.boundary_extension_right += 1;
    }
}

fn try_extend_repeat_core_pair(
    pair: &mut RepeatCorePair,
    record: &FastaRecord,
    side: &str,
) -> bool {
    let reference_len = record.sequence.len();
    let Some(pos_a) = repeat_copy_extension_position(&pair.copies[0], side, reference_len) else {
        return false;
    };
    let Some(pos_b) = repeat_copy_extension_position(&pair.copies[1], side, reference_len) else {
        return false;
    };
    let base_a = oriented_base_at(record, pos_a, pair.copies[0].strand);
    let base_b = oriented_base_at(record, pos_b, pair.copies[1].strand);
    if !base_a.eq_ignore_ascii_case(&base_b) {
        return false;
    }
    let copy_a = extended_repeat_copy(&pair.copies[0], pos_a, side);
    let copy_b = extended_repeat_copy(&pair.copies[1], pos_b, side);
    if interval_overlap(copy_a.start, copy_a.end, copy_b.start, copy_b.end) > 0 {
        return false;
    }
    pair.copies = vec![copy_a, copy_b];
    pair.length += 1;
    true
}

fn repeat_copy_extension_position(
    copy: &RepeatCopy,
    side: &str,
    reference_len: usize,
) -> Option<usize> {
    if side == "left" {
        if copy.strand == '+' {
            (copy.start > 1).then_some(copy.start - 1)
        } else {
            (copy.end < reference_len).then_some(copy.end + 1)
        }
    } else if copy.strand == '+' {
        (copy.end < reference_len).then_some(copy.end + 1)
    } else {
        (copy.start > 1).then_some(copy.start - 1)
    }
}

fn extended_repeat_copy(copy: &RepeatCopy, position: usize, side: &str) -> RepeatCopy {
    let mut updated = copy.clone();
    if side == "left" {
        if updated.strand == '+' {
            updated.start = position;
        } else {
            updated.end = position;
        }
    } else if updated.strand == '+' {
        updated.end = position;
    } else {
        updated.start = position;
    }
    updated
}

fn oriented_base_at(record: &FastaRecord, position: usize, strand: char) -> char {
    let base = record
        .sequence
        .as_bytes()
        .get(position.saturating_sub(1))
        .copied()
        .unwrap_or(b'N') as char;
    if strand == '-' {
        complement_base(base)
    } else {
        base
    }
}

fn assign_repeat_nodes_to_core_pairs(
    merged: &Gfa,
    node_classes: &BTreeMap<String, String>,
    mapping: &Mapping,
    core_pairs: &[RepeatCorePair],
) -> BTreeMap<String, RepeatCorePair> {
    let mut assignments = BTreeMap::new();
    let mut used = HashSet::new();
    for name in &merged.order {
        if node_classes[name] != "repeat_node" {
            continue;
        }
        let node_mapping = &mapping.by_node[name];
        let relaxed_anchor = node_mapping.accepted_hits.len() < 2;
        let mut candidates = Vec::new();
        for (index, pair) in core_pairs.iter().enumerate() {
            if used.contains(&index) {
                continue;
            }
            let Some(score) = repeat_pair_mapping_score(node_mapping, pair, mapping.reference_len)
            else {
                continue;
            };
            candidates.push((index, score));
        }
        let selected = if relaxed_anchor {
            // A relaxed repeat anchor must identify exactly one self-BLAST core pair.
            // Ambiguous pairs remain fail-closed and use the existing fallback plan.
            (candidates.len() == 1).then(|| candidates[0].0)
        } else {
            candidates
                .into_iter()
                .max_by_key(|(_, score)| *score)
                .map(|(index, _)| index)
        };
        if let Some(index) = selected {
            used.insert(index);
            assignments.insert(name.clone(), core_pairs[index].clone());
        }
    }
    assignments
}

fn repeat_assignment_hits(mapping: &NodeMapping, reference_len: usize) -> Vec<PafHit> {
    if mapping.accepted_hits.len() >= 2 {
        return mapping.accepted_hits.clone();
    }
    dedupe_hits(
        mapping
            .all_hits
            .iter()
            .filter(|hit| {
                hit.aligned_fraction() >= REPEAT_ASSIGNMENT_MIN_ALIGNED_FRACTION
                    && hit.substitution_identity() >= REPEAT_ASSIGNMENT_MIN_SUBSTITUTION_IDENTITY
            })
            .cloned()
            .collect(),
        reference_len,
    )
}

fn repeat_pair_mapping_score(
    mapping: &NodeMapping,
    pair: &RepeatCorePair,
    reference_len: usize,
) -> Option<usize> {
    if pair.copies.len() != 2 {
        return None;
    }
    let assignment_hits = repeat_assignment_hits(mapping, reference_len);
    if assignment_hits.len() < 2 {
        return None;
    }
    let required_overlap = pair.length.saturating_mul(3).div_ceil(4);
    let mut best = None;
    for first_index in 0..assignment_hits.len() {
        for second_index in 0..assignment_hits.len() {
            if first_index == second_index {
                continue;
            }
            let first_overlap = circular_hit_interval_overlap(
                &assignment_hits[first_index],
                &pair.copies[0],
                reference_len,
            );
            let second_overlap = circular_hit_interval_overlap(
                &assignment_hits[second_index],
                &pair.copies[1],
                reference_len,
            );
            if first_overlap < required_overlap || second_overlap < required_overlap {
                continue;
            }
            let score = first_overlap + second_overlap;
            if best.is_none_or(|current| score > current) {
                best = Some(score);
            }
        }
    }
    best
}

fn build_repeat_core_gaps(
    repeat_assignments: &BTreeMap<String, RepeatCorePair>,
    reference_len: usize,
) -> Vec<RepeatGap> {
    let mut copies = Vec::new();
    for (node, pair) in repeat_assignments {
        for (copy_index, copy) in pair.copies.iter().enumerate() {
            copies.push((node.clone(), copy_index + 1, copy.clone()));
        }
    }
    copies.sort_by_key(|(_, _, copy)| copy.start);
    let mut gaps = Vec::new();
    if copies.is_empty() {
        return gaps;
    }
    for index in 0..copies.len() {
        let (left_node, left_copy_index, left_copy) = &copies[index];
        let (right_node, right_copy_index, right_copy) = &copies[(index + 1) % copies.len()];
        let parts = repeat_gap_parts(left_copy, right_copy, reference_len);
        gaps.push(RepeatGap {
            left_repeat_node: left_node.clone(),
            left_copy_index: *left_copy_index,
            right_repeat_node: right_node.clone(),
            right_copy_index: *right_copy_index,
            parts,
        });
    }
    gaps
}

fn repeat_gap_parts(
    left_copy: &RepeatCopy,
    right_copy: &RepeatCopy,
    reference_len: usize,
) -> Vec<(usize, usize)> {
    let start = left_copy.end + 1;
    let end = right_copy.start.saturating_sub(1);
    if start <= end {
        vec![(start, end)]
    } else {
        let mut parts = Vec::new();
        if start <= reference_len {
            parts.push((start, reference_len));
        }
        if end >= 1 {
            parts.push((1, end));
        }
        parts
    }
}

fn assign_single_copy_nodes_to_gaps(
    merged: &Gfa,
    node_classes: &BTreeMap<String, String>,
    mapping: &Mapping,
    gaps: &[RepeatGap],
) -> BTreeMap<String, RepeatGap> {
    let mut assignments = BTreeMap::new();
    let mut used = HashSet::new();
    for name in &merged.order {
        if node_classes[name] != "single_copy_node" {
            continue;
        }
        let mut best_index = None;
        let mut best_score = 0usize;
        for (index, gap) in gaps.iter().enumerate() {
            if used.contains(&index) {
                continue;
            }
            let score: usize = mapping.by_node[name]
                .all_hits
                .iter()
                .map(|hit| hit_overlap_with_parts(hit, &gap.parts))
                .sum();
            if best_index.is_none() || score > best_score {
                best_index = Some(index);
                best_score = score;
            }
        }
        if let Some(index) = best_index.filter(|_| best_score > 0) {
            used.insert(index);
            assignments.insert(name.clone(), gaps[index].clone());
        }
    }
    assignments
}

fn choose_repeat_copy_for_node(
    name: &str,
    pair: &RepeatCorePair,
    mapping: &Mapping,
) -> (RepeatCopy, char) {
    let mut best_copy = pair.copies[0].clone();
    let mut best_strand = pair.copies[0].strand;
    let mut best_score = 0usize;
    for hit in &mapping.by_node[name].all_hits {
        for copy in &pair.copies {
            let score = hit_interval_overlap(hit, copy);
            if score > best_score {
                best_score = score;
                best_copy = copy.clone();
                best_strand = hit.strand;
            }
        }
    }
    (best_copy, best_strand)
}

fn best_orientation_for_parts(name: &str, parts: &[(usize, usize)], mapping: &Mapping) -> char {
    let mut best_strand = '+';
    let mut best_score = 0usize;
    for hit in &mapping.by_node[name].all_hits {
        let score = hit_overlap_with_parts(hit, parts);
        if score > best_score {
            best_score = score;
            best_strand = hit.strand;
        }
    }
    best_strand
}

fn hit_interval_overlap(hit: &PafHit, copy: &RepeatCopy) -> usize {
    interval_overlap(hit.target_start + 1, hit.target_end, copy.start, copy.end)
}

fn circular_hit_interval_overlap(hit: &PafHit, copy: &RepeatCopy, reference_len: usize) -> usize {
    folded_interval_parts(&hit.folded(reference_len), reference_len)
        .into_iter()
        .map(|(start, end)| interval_overlap(start, end, copy.start, copy.end))
        .sum()
}

fn hit_overlap_with_parts(hit: &PafHit, parts: &[(usize, usize)]) -> usize {
    parts
        .iter()
        .map(|(start, end)| interval_overlap(hit.target_start + 1, hit.target_end, *start, *end))
        .sum()
}

fn synthetic_hits_for_repeat_copies(
    name: &str,
    node_sequence: &str,
    record: &FastaRecord,
    copies: &[RepeatCopy],
) -> Vec<PafHit> {
    let mut hits = Vec::new();
    for copy in copies {
        let plus_sequence = record_interval_sequence(record, copy.start, copy.end, '+');
        let minus_sequence = reverse_complement(&plus_sequence);
        let strand = if node_sequence == plus_sequence {
            '+'
        } else if node_sequence == minus_sequence {
            '-'
        } else {
            copy.strand
        };
        hits.push(make_synthetic_hit(
            name,
            node_sequence.len(),
            1,
            node_sequence.len(),
            copy.start,
            copy.end,
            strand,
        ));
    }
    hits
}

fn synthetic_hits_for_parts(
    name: &str,
    node_len: usize,
    record: &FastaRecord,
    parts: &[(usize, usize)],
    strand: char,
) -> Vec<PafHit> {
    let mut hits = Vec::new();
    let mut query_pos = 1usize;
    let mut ordered_parts = parts.to_vec();
    if strand == '-' {
        ordered_parts.reverse();
    }
    for (start, end) in ordered_parts {
        let part_len = end - start + 1;
        hits.push(make_synthetic_hit(
            name,
            node_len,
            query_pos,
            query_pos + part_len - 1,
            start,
            end,
            strand,
        ));
        query_pos += part_len;
    }
    // Keep reference id in use so future refactoring does not accidentally change semantics.
    let _ = &record.id;
    hits
}

fn make_synthetic_hit(
    name: &str,
    node_len: usize,
    query_start: usize,
    query_end: usize,
    target_start: usize,
    target_end: usize,
    strand: char,
) -> PafHit {
    PafHit {
        query_name: name.to_string(),
        query_length: node_len,
        query_start: query_start - 1,
        query_end,
        strand,
        target_start: target_start - 1,
        target_end,
        matches: query_end - query_start + 1,
        block_length: query_end - query_start + 1,
        mapq: 60,
        substitution_matches: query_end - query_start + 1,
        substitution_mismatches: 0,
        primary: true,
        exact_flanks: false,
    }
}

fn verified_sequence_from_parts(
    record: &FastaRecord,
    parts: &[(usize, usize)],
    strand: char,
) -> String {
    let mut sequence = String::new();
    for (start, end) in parts {
        sequence.push_str(&record_interval_sequence(record, *start, *end, '+'));
    }
    if strand == '-' {
        reverse_complement(&sequence)
    } else {
        sequence
    }
}

fn extract_verified_sequence(record: &FastaRecord, hit: &PafHit) -> String {
    let mut sequence = extract_circular(&record.sequence, hit.target_start, hit.target_end);
    if hit.strand == '-' {
        sequence = reverse_complement(&sequence);
    }
    sequence
}

fn record_interval_sequence(
    record: &FastaRecord,
    start: usize,
    end: usize,
    strand: char,
) -> String {
    let sequence = record
        .sequence
        .chars()
        .skip(start.saturating_sub(1))
        .take(end.saturating_sub(start) + 1)
        .collect::<String>();
    if strand == '-' {
        reverse_complement(&sequence)
    } else {
        sequence
    }
}

fn interval_overlap(start_a: usize, end_a: usize, start_b: usize, end_b: usize) -> usize {
    let start = start_a.max(start_b);
    let end = end_a.min(end_b);
    if end < start {
        0
    } else {
        end - start + 1
    }
}

fn sorted_pair(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

fn write_self_blast_repeat_core_debug(
    path: &Path,
    pairs: &[RepeatCorePair],
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "pair_index\tlength\tsource_blast_length\tsource_pident\tsource_mismatch\tsource_gapopen\tcopy_index\tstart\tend\tstrand\tboundary_extension_left\tboundary_extension_right")?;
    for (pair_index, pair) in pairs.iter().enumerate() {
        for (copy_index, copy) in pair.copies.iter().enumerate() {
            writeln!(
                out,
                "{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                pair_index + 1,
                pair.length,
                pair.source_blast_length,
                pair.source_pident,
                pair.source_mismatch,
                pair.source_gapopen,
                copy_index + 1,
                copy.start,
                copy.end,
                copy.strand,
                pair.boundary_extension_left,
                pair.boundary_extension_right
            )?;
        }
    }
    Ok(())
}

fn write_rebuild_extract(
    path: &Path,
    subgraph: &str,
    gfa: &Gfa,
    mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    source_rows: &[VerifiedNodeSourceRow],
    coverage_rows: &[CoverageRow],
) -> Result<(), OrgraftError> {
    let degrees = gfa.degrees();
    let source_by_node: BTreeMap<String, VerifiedNodeSourceRow> = source_rows
        .iter()
        .map(|row| (row.node.clone(), row.clone()))
        .collect();
    let coverage_by_node: BTreeMap<String, CoverageRow> = coverage_rows
        .iter()
        .map(|row| (row.node.clone(), row.clone()))
        .collect();
    let mut rows: Vec<ExtractRow> = Vec::new();
    for name in &gfa.order {
        let segment = &gfa.segments[name];
        let Some(node) = mapping.by_node.get(name) else {
            rows.push(ExtractRow::unmapped(
                subgraph,
                name,
                segment.sequence.len(),
                node_classes[name].as_str(),
                degrees.get(name).copied().unwrap_or(0),
            ));
            continue;
        };
        for (copy_index, hit) in node.accepted_hits.iter().enumerate() {
            let interval = hit.folded(mapping.reference_len);
            let source = source_by_node.get(name);
            let coverage = coverage_by_node.get(name);
            let (node_start, node_end) = if hit.strand == '-' {
                (hit.query_end, hit.query_start + 1)
            } else {
                (hit.query_start + 1, hit.query_end)
            };
            let mut push_interval = |linear_start: usize, linear_end: usize, notes: &str| {
                rows.push(ExtractRow {
                    subgraph: subgraph.to_string(),
                    row_type: "node_copy".to_string(),
                    linear_start: linear_start.to_string(),
                    linear_end: linear_end.to_string(),
                    node: name.clone(),
                    node_start,
                    node_end,
                    strand: hit.strand,
                    node_length: segment.sequence.len(),
                    node_class: node_classes[name].clone(),
                    copy_index: copy_index + 1,
                    source_kind: source
                        .map(|row| row.source.clone())
                        .unwrap_or_else(|| "polished_fasta_minimap2".to_string()),
                    source_node_class: source
                        .map(|row| row.node_class.clone())
                        .unwrap_or_else(|| ".".to_string()),
                    source_length: source
                        .map(|row| row.length.to_string())
                        .unwrap_or_else(|| ".".to_string()),
                    source_parts: source
                        .map(|row| row.parts.clone())
                        .unwrap_or_else(|| ".".to_string()),
                    source_notes: source
                        .map(|row| row.notes.clone())
                        .unwrap_or_else(|| ".".to_string()),
                    depth_mean: optional_f64(coverage.and_then(|row| row.depth_mean)),
                    coverage_intervals: coverage
                        .map(|row| row.intervals.clone())
                        .unwrap_or_else(|| ".".to_string()),
                    coverage_source: coverage
                        .map(|row| row.coverage_source.clone())
                        .unwrap_or_else(|| ".".to_string()),
                    degree: degrees.get(name).copied().unwrap_or(0),
                    identity: format!("{:.6}", hit.identity()),
                    aligned_fraction: format!("{:.6}", hit.aligned_fraction()),
                    mapq: hit.mapq.to_string(),
                    notes: notes.to_string(),
                });
            };
            for (part_index, (start, end)) in interval.parts.iter().enumerate() {
                let notes =
                    if hit.target_end.saturating_sub(hit.target_start) > mapping.reference_len {
                        format!(
                            "wrapped_circular_interval;part={}/{}",
                            part_index + 1,
                            interval.parts.len()
                        )
                    } else if interval.parts.len() > 1 {
                        "wrapped_circular_interval".to_string()
                    } else {
                        ".".to_string()
                    };
                push_interval(*start, *end, &notes);
            }
        }
        if node.accepted_hits.is_empty() {
            rows.push(ExtractRow::unmapped(
                subgraph,
                name,
                segment.sequence.len(),
                node_classes[name].as_str(),
                degrees.get(name).copied().unwrap_or(0),
            ));
        }
    }
    rows.sort_by(|a, b| {
        parse_sortable_usize(&a.linear_start)
            .cmp(&parse_sortable_usize(&b.linear_start))
            .then_with(|| {
                parse_sortable_usize(&a.linear_end).cmp(&parse_sortable_usize(&b.linear_end))
            })
            .then_with(|| natural_cmp(&a.node, &b.node))
            .then_with(|| a.copy_index.cmp(&b.copy_index))
    });
    let mut out = File::create(path)?;
    writeln!(out, "subgraph\trow_type\tlinear_start\tlinear_end\tnode\tnode_start\tnode_end\tstrand\tnode_length\tnode_class\tcopy_index\tsource_kind\tsource_node_class\tsource_length\tsource_parts\tsource_notes\tdepth_mean\tcoverage_intervals\tcoverage_source\tdegree\tidentity\taligned_fraction\tmapq\tnotes")?;
    for row in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.subgraph,
            row.row_type,
            row.linear_start,
            row.linear_end,
            row.node,
            row.node_start,
            row.node_end,
            row.strand,
            row.node_length,
            row.node_class,
            row.copy_index,
            row.source_kind,
            row.source_node_class,
            row.source_length,
            row.source_parts,
            row.source_notes,
            row.depth_mean,
            row.coverage_intervals,
            row.coverage_source,
            row.degree,
            row.identity,
            row.aligned_fraction,
            row.mapq,
            row.notes
        )?;
    }
    Ok(())
}

#[derive(Debug)]
struct ExtractRow {
    subgraph: String,
    row_type: String,
    linear_start: String,
    linear_end: String,
    node: String,
    node_start: usize,
    node_end: usize,
    strand: char,
    node_length: usize,
    node_class: String,
    copy_index: usize,
    source_kind: String,
    source_node_class: String,
    source_length: String,
    source_parts: String,
    source_notes: String,
    depth_mean: String,
    coverage_intervals: String,
    coverage_source: String,
    degree: usize,
    identity: String,
    aligned_fraction: String,
    mapq: String,
    notes: String,
}

impl ExtractRow {
    fn unmapped(
        subgraph: &str,
        node: &str,
        node_length: usize,
        node_class: &str,
        degree: usize,
    ) -> Self {
        Self {
            subgraph: subgraph.to_string(),
            row_type: "node_unmapped".to_string(),
            linear_start: ".".to_string(),
            linear_end: ".".to_string(),
            node: node.to_string(),
            node_start: 0,
            node_end: 0,
            strand: '.',
            node_length,
            node_class: node_class.to_string(),
            copy_index: 0,
            source_kind: ".".to_string(),
            source_node_class: ".".to_string(),
            source_length: ".".to_string(),
            source_parts: ".".to_string(),
            source_notes: "no_alignment".to_string(),
            depth_mean: ".".to_string(),
            coverage_intervals: ".".to_string(),
            coverage_source: ".".to_string(),
            degree,
            identity: ".".to_string(),
            aligned_fraction: ".".to_string(),
            mapq: ".".to_string(),
            notes: "no_alignment".to_string(),
        }
    }
}

fn parse_sortable_usize(value: &str) -> usize {
    value.parse::<usize>().unwrap_or(usize::MAX)
}

fn export_gfa_reference_images(
    paths: &OutputPaths,
    image_reference_fasta: &Path,
    gfa_editor_cli: &Path,
    soft_paths: &HashMap<String, PathBuf>,
) -> Vec<GfaImageExport> {
    let mut rows = Vec::new();
    for (format, output_path) in [
        ("pdf", paths.verified_pdf.as_path()),
        ("svg", paths.verified_svg.as_path()),
    ] {
        rows.push(run_gfa_editor_image(
            gfa_editor_cli,
            soft_paths,
            &paths.verified_gfa,
            output_path,
            image_reference_fasta,
            format,
        ));
    }
    rows
}

fn skipped_gfa_reference_images(paths: &OutputPaths, reason: &str) -> Vec<GfaImageExport> {
    [
        ("pdf", paths.verified_pdf.as_path()),
        ("svg", paths.verified_svg.as_path()),
    ]
    .into_iter()
    .map(|(format, output_path)| GfaImageExport {
        format: format.to_string(),
        output: output_path.to_path_buf(),
        command: ".".to_string(),
        status: "skipped_missing_gfa_editor_cli".to_string(),
        stdout: String::new(),
        stderr: reason.to_string(),
    })
    .collect()
}

fn write_run_report(
    path: &Path,
    out_dir: &Path,
    paths: &OutputPaths,
    options: &RebuildOptions,
    edited_gfa: &Path,
    polished_fasta: &Path,
    reads: Option<&Path>,
    image_reference_fasta: Option<&Path>,
    validate_data_dir: Option<&Path>,
    merge_mode: &str,
    repeat_resolution: Option<&RepeatResolutionOutcome>,
    sequence_projection_method: &str,
    minimap2: &Path,
    blastn: &Path,
    gfa_editor_cli: Option<&Path>,
    mapping: &Mapping,
    image_exports: &[GfaImageExport],
    started: &Instant,
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "section\tkey\tvalue")?;
    let proof=paths.logs_dir.join(format!("rebuild_{}_ordinary_source.tsv",options.subgraph));
    if let Some(source) = options.ordinary_source_dir.as_ref().filter(|_| proof.is_file()) {
        writeln!(out, "input\tordinary_source_dir\t{}", fs::canonicalize(source)?.display())?;
        writeln!(out, "output\tordinary_source_proof\t{}", fs::canonicalize(&proof)?.display())?;
        writeln!(out, "output\tordinary_source_proof_sha256\t{}", sha256_hex(&fs::read(proof)?))?;
    }
    let mut rows = vec![
        ("run", "status", "completed".to_string()),
        ("run", "organelle", options.organelle.clone()),
        ("run", "subgraph", options.subgraph.clone()),
        ("run", "threads", options.threads.to_string()),
        (
            "run",
            "runtime_seconds",
            format!("{:.6}", started.elapsed().as_secs_f64()),
        ),
        ("input", "edited_gfa", edited_gfa.display().to_string()),
        (
            "input",
            "structure_corrections",
            options
                .structure_corrections
                .as_ref()
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "input",
            "checked_gfa",
            options
                .checked_gfa
                .as_ref()
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "input",
            "resolve_details",
            options
                .resolve_details
                .as_ref()
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "input",
            "polished_fasta",
            polished_fasta.display().to_string(),
        ),
        (
            "input",
            "reads",
            reads
                .map(|path| path.display().to_string())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "input",
            "image_reference_fasta",
            image_reference_fasta
                .map(|path| path.display().to_string())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "input",
            "validate_data_dir",
            validate_data_dir
                .map(|path| path.display().to_string())
                .unwrap_or_else(|| ".".to_string()),
        ),
        (
            "output",
            "rebuild_gfa",
            relative_output_path(out_dir, &paths.verified_gfa),
        ),
        (
            "output",
            "rebuild_fasta",
            relative_output_path(out_dir, &paths.verified_fasta),
        ),
        (
            "output",
            "rebuild_nodes_fasta",
            relative_output_path(out_dir, &paths.verified_nodes_fasta),
        ),
        (
            "output",
            "extract_table",
            relative_output_path(out_dir, &paths.extract_report),
        ),
        (
            "output",
            "run_report",
            relative_output_path(out_dir, &paths.run_report),
        ),
        (
            "output",
            "result_stats",
            relative_output_path(out_dir, &paths.result_stats),
        ),
        (
            "output",
            "repeat_path_support",
            relative_output_path(out_dir, &paths.repeat_path_support),
        ),
        (
            "output",
            "repeat_resolution",
            relative_output_path(out_dir, &paths.repeat_resolution),
        ),
        (
            "file_description",
            "rebuild_fasta",
            "complete verified/polished linear sequence for the subgraph".to_string(),
        ),
        (
            "file_description",
            "rebuild_nodes_fasta",
            "FASTA records extracted from each S node in the rebuilt GFA".to_string(),
        ),
        ("method", "merge_mode", merge_mode.to_string()),
        (
            "method",
            "repeat_resolution",
            repeat_resolution
                .map(|outcome| {
                    format!(
                        "spanning_read_constrained:{}:source={}:local={}",
                        outcome.decision.repeat_node,
                        outcome.selected_candidate,
                        outcome.local_candidate,
                    )
                })
                .unwrap_or_else(|| "disabled".to_string()),
        ),
        (
            "method",
            "sequence_projection",
            sequence_projection_method.to_string(),
        ),
        ("tool", "minimap2", minimap2.display().to_string()),
        ("tool", "blastn", blastn.display().to_string()),
        (
            "tool",
            "gfa_editor_cli",
            gfa_editor_cli
                .map(|path| path.display().to_string())
                .unwrap_or_else(|| {
                    if image_reference_fasta.is_some() {
                        "skipped_missing_gfa_editor_cli".to_string()
                    } else {
                        "skipped_no_image_reference".to_string()
                    }
                }),
        ),
        ("tool", "minimap2_command", mapping.command.clone()),
        (
            "tool",
            "minimap2_stderr",
            mapping.stderr.replace('\n', "\\n"),
        ),
    ];
    if let (Some(trial), Some(evidence)) = (&options.joint_ir_trial_gfa, &options.joint_ir_evidence) {
        for (key, hash_key, source) in [
            ("joint_ir_trial_gfa", "joint_ir_trial_gfa_sha256", trial),
            ("joint_ir_evidence", "joint_ir_evidence_sha256", evidence),
        ] {
            rows.push(("input", key, fs::canonicalize(source)?.display().to_string()));
            rows.push(("input", hash_key, file_sha(source)?));
        }
    }
    for row in image_exports.iter().filter(|row| row.status == "written") {
        let key = match row.format.as_str() {
            "pdf" => "rebuild_pdf",
            "svg" => "rebuild_svg",
            other => other,
        };
        rows.push(("output", key, relative_output_path(out_dir, &row.output)));
    }
    for (section, key, value) in rows {
        writeln!(out, "{section}\t{key}\t{}", value.replace('\n', "\\n"))?;
    }
    for row in image_exports {
        writeln!(
            out,
            "image\t{}_output\t{}",
            row.format,
            relative_output_path(out_dir, &row.output)
        )?;
        writeln!(out, "image\t{}_status\t{}", row.format, row.status)?;
        writeln!(
            out,
            "image\t{}_command\t{}",
            row.format,
            row.command.replace('\n', "\\n")
        )?;
        writeln!(
            out,
            "image\t{}_stdout\t{}",
            row.format,
            row.stdout.replace('\n', "\\n")
        )?;
        writeln!(
            out,
            "image\t{}_stderr\t{}",
            row.format,
            row.stderr.replace('\n', "\\n")
        )?;
    }
    Ok(())
}

fn write_result_stats(
    path: &Path,
    subgraph: &str,
    summary: &[(&str, String)],
    raw_stats: &GraphStats,
    merged_stats: &GraphStats,
    verified_stats: &GraphStats,
    consistency: &CoordinateConsistency,
    coverage_rows: &[CoverageRow],
    repeat_rows: &[RepeatPathSupportRow],
) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    writeln!(out, "section\tsubgraph\titem\tmetric\tvalue\textra")?;
    for (key, value) in summary {
        writeln!(
            out,
            "summary\t{subgraph}\t.\t{key}\t{}\t.",
            value.replace('\n', "\\n")
        )?;
    }
    for (name, stats) in [
        ("raw", raw_stats),
        ("merged", merged_stats),
        ("verified", verified_stats),
    ] {
        writeln!(
            out,
            "graph\t{subgraph}\t{name}\tsegments\t{}\t.",
            stats.segments
        )?;
        writeln!(out, "graph\t{subgraph}\t{name}\tlinks\t{}\t.", stats.links)?;
        writeln!(out, "graph\t{subgraph}\t{name}\tbp\t{}\t.", stats.bp)?;
    }
    writeln!(
        out,
        "consistency\t{subgraph}\tpolished\tlength\t{}\thash64={}",
        consistency.polished_length, consistency.polished_hash64
    )?;
    writeln!(
        out,
        "consistency\t{subgraph}\tsegments\tmapped\t{}\tunmapped={}",
        consistency.mapped_segments, consistency.unmapped_segments
    )?;
    writeln!(
        out,
        "consistency\t{subgraph}\tcoverage\tcovered_bases\t{}\tmulti_covered_bases={};gap_bases={};coverage_fraction={:.6};linear_tiling={}",
        consistency.covered_bases,
        consistency.multi_covered_bases,
        consistency.gap_bases,
        consistency.coverage_fraction,
        consistency.linear_tiling_status()
    )?;
    writeln!(
        out,
        "consistency\t{subgraph}\tgraph_segments\tverified\t{}\traw={};merged={}",
        consistency.verified_segments, consistency.raw_segments, consistency.merged_segments
    )?;
    for row in coverage_rows {
        writeln!(
            out,
            "node_depth\t{subgraph}\t{}\tdepth_mean\t{}\tnode_class={};covered_bases={};intervals={};source={};coverage_path={}",
            row.node,
            optional_f64(row.depth_mean),
            row.node_class,
            optional_usize(row.covered_bases),
            row.intervals,
            row.coverage_source,
            row.coverage_path
        )?;
    }
    for row in repeat_rows {
        writeln!(
            out,
            "repeat_path\t{subgraph}\t{}\tpath_support\t{}\tstatus={};method={};left={};right={};ratio={};left_edge={};right_edge={};left_support={};right_support={}",
            row.repeat_node,
            optional_usize(row.path_support),
            row.repeat_status,
            row.support_method,
            row.left_endpoint,
            row.right_endpoint,
            optional_f64(row.path_ratio),
            optional_usize(row.left_edge_index),
            optional_usize(row.right_edge_index),
            optional_f64(row.left_support),
            optional_f64(row.right_support)
        )?;
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct CoverageRow {
    node: String,
    node_class: String,
    coverage_source: String,
    depth_mean: Option<f64>,
    depth_bases: Option<f64>,
    covered_bases: Option<usize>,
    interval_count: usize,
    intervals: String,
    coverage_path: String,
}

fn compute_node_remapped_coverage(
    mapping: &Mapping,
    node_classes: &BTreeMap<String, String>,
    validate_data_dir: Option<&Path>,
) -> Result<Vec<CoverageRow>, OrgraftError> {
    let coverage_path = validate_data_dir
        .map(|dir| dir.join("sv_coverage.tsv"))
        .unwrap_or_default();
    let depths = read_remapped_depths(&coverage_path)?;
    let mut rows = Vec::new();
    for (name, node) in &mapping.by_node {
        let intervals: Vec<(usize, usize)> = node
            .accepted_hits
            .iter()
            .flat_map(|hit| hit.folded(mapping.reference_len).parts)
            .collect();
        // Depths refer to physical positions in the one-copy authoritative FASTA.
        // Fold complete hit spans before taking the union, so a crossing hit keeps
        // its leading positions and repeated turns do not multiply their depth.
        let merged = merge_intervals(&intervals);
        if let Some((depth_mean, depth_bases, covered_bases)) =
            coverage_for_intervals(&depths, &merged)
        {
            rows.push(CoverageRow {
                node: name.clone(),
                node_class: node_classes[name].clone(),
                coverage_source: "remapped_FL_to_verified_fasta".to_string(),
                depth_mean: Some(depth_mean),
                depth_bases: Some(depth_bases),
                covered_bases: Some(covered_bases),
                interval_count: merged.len(),
                intervals: intervals_text(&merged),
                coverage_path: coverage_path.display().to_string(),
            });
        } else {
            rows.push(CoverageRow {
                node: name.clone(),
                node_class: node_classes[name].clone(),
                coverage_source: "remapped_FL_to_verified_fasta_unavailable".to_string(),
                depth_mean: None,
                depth_bases: None,
                covered_bases: None,
                interval_count: 0,
                intervals: ".".to_string(),
                coverage_path: coverage_path.display().to_string(),
            });
        }
    }
    Ok(rows)
}

fn read_remapped_depths(path: &Path) -> Result<Vec<f64>, OrgraftError> {
    if !path.exists() {
        return Ok(Vec::new());
    }
    let text = fs::read_to_string(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
    })?;
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return Ok(Vec::new());
    };
    let columns: Vec<&str> = header.split('\t').collect();
    let depth_index = columns
        .iter()
        .position(|column| *column == "fl_depth")
        .or_else(|| columns.iter().position(|column| *column == "depth"))
        .unwrap_or_else(|| columns.len().saturating_sub(1));
    let mut depths = vec![0.0];
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        let value = fields
            .get(depth_index)
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(0.0);
        depths.push(value);
    }
    Ok(depths)
}

fn merge_intervals(intervals: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut sorted: Vec<(usize, usize)> = intervals
        .iter()
        .copied()
        .filter(|(start, end)| start <= end)
        .collect();
    sorted.sort();
    let mut merged: Vec<(usize, usize)> = Vec::new();
    for (start, end) in sorted {
        if let Some(last) = merged.last_mut() {
            if start <= last.1 + 1 {
                last.1 = last.1.max(end);
                continue;
            }
        }
        merged.push((start, end));
    }
    merged
}

fn coverage_for_intervals(
    depths: &[f64],
    intervals: &[(usize, usize)],
) -> Option<(f64, f64, usize)> {
    if depths.is_empty() || intervals.is_empty() {
        return None;
    }
    let max_position = depths.len().saturating_sub(1);
    let mut depth_bases = 0.0;
    let mut covered_bases = 0usize;
    for (start, end) in intervals {
        let start = (*start).max(1);
        let end = (*end).min(max_position);
        if end < start {
            continue;
        }
        for value in &depths[start..=end] {
            depth_bases += *value;
        }
        covered_bases += end - start + 1;
    }
    if covered_bases == 0 {
        None
    } else {
        Some((
            depth_bases / covered_bases as f64,
            depth_bases,
            covered_bases,
        ))
    }
}

fn intervals_text(intervals: &[(usize, usize)]) -> String {
    if intervals.is_empty() {
        ".".to_string()
    } else {
        intervals
            .iter()
            .map(|(start, end)| format!("{start}-{end}"))
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn annotate_gfa_with_remapped_coverage(gfa: &mut Gfa, rows: &[CoverageRow]) {
    let by_node: BTreeMap<&str, &CoverageRow> =
        rows.iter().map(|row| (row.node.as_str(), row)).collect();
    for name in &gfa.order {
        let Some(row) = by_node.get(name.as_str()) else {
            continue;
        };
        let Some(segment) = gfa.segments.get_mut(name) else {
            continue;
        };
        let mut tags = segment.tags.clone();
        if let (Some(depth_mean), Some(depth_bases)) = (row.depth_mean, row.depth_bases) {
            tags = replace_tags(
                &tags,
                &[
                    format!("DP:f:{depth_mean:.6}"),
                    "CM:Z:remapped_FL_to_verified_fasta".to_string(),
                    format!("RI:i:{}", row.interval_count),
                    format!("RB:f:{depth_bases:.3}"),
                ],
            );
        } else {
            tags = replace_tag(
                &tags,
                "CM:Z:remapped_FL_to_verified_fasta_unavailable".to_string(),
            );
        }
        segment.tags = tags;
    }
}

#[derive(Debug, Clone)]
pub(crate) struct RepeatPathSupportRow {
    pub(crate) repeat_node: String,
    pub(crate) repeat_status: String,
    pub(crate) support_method: String,
    pub(crate) left_edge_index: Option<usize>,
    pub(crate) right_edge_index: Option<usize>,
    pub(crate) left_endpoint: String,
    pub(crate) right_endpoint: String,
    pub(crate) left_physical_endpoint: String,
    pub(crate) right_physical_endpoint: String,
    pub(crate) left_support: Option<f64>,
    pub(crate) right_support: Option<f64>,
    pub(crate) path_support: Option<usize>,
    pub(crate) path_ratio: Option<f64>,
    pub(crate) supporting_read_ids: Vec<String>,
}

#[derive(Debug, Clone)]
struct RepeatPathEdge {
    index: usize,
    endpoint_name: String,
    endpoint_orient: char,
    endpoint_side: char,
    support: f64,
    overlap: usize,
}

impl RepeatPathEdge {
    fn label(&self) -> String {
        oriented_node(&self.endpoint_name, self.endpoint_orient)
    }

    fn physical_label(&self) -> String {
        oriented_node(&self.endpoint_name, self.endpoint_side)
    }
}

#[derive(Debug, Clone)]
struct RepeatPathTemplate {
    id: String,
    repeat_node: String,
    row_index: usize,
    sequence: String,
    required_start: usize,
    required_end: usize,
}

#[derive(Debug, Clone)]
struct RepeatPathAlignment {
    read_id: String,
    target_id: String,
    query_start: usize,
    query_end: usize,
    target_start: usize,
    target_end: usize,
    matches: usize,
    mapq: usize,
    primary: bool,
}

impl RepeatPathAlignment {
    fn identity(&self) -> f64 {
        let aligned_span = self.query_span().min(self.target_span());
        if aligned_span == 0 {
            0.0
        } else {
            self.matches as f64 / aligned_span as f64
        }
    }

    fn query_span(&self) -> usize {
        self.query_end.saturating_sub(self.query_start)
    }

    fn target_span(&self) -> usize {
        self.target_end.saturating_sub(self.target_start)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct RepeatResolutionDecision {
    pub(crate) repeat_node: String,
    pub(crate) pairing_signature: String,
    pub(crate) selected_paths: Vec<String>,
    pub(crate) selected_support: usize,
    pub(crate) alternative_support: usize,
}

#[derive(Debug, Clone)]
struct RepeatResolutionOutcome {
    decision: RepeatResolutionDecision,
    // Stable ID from the first constrained rebuild on the original candidate set.
    selected_candidate: String,
    // ID after re-enumerating candidates against this rebuild's corrected sequence.
    local_candidate: String,
    candidate_count: usize,
    filtered_candidate_count: usize,
    selection_method: String,
    score: Option<f64>,
    continuous_fraction: Option<f64>,
    selected_order: String,
    ordered_nodes: Vec<String>,
    oriented_nodes: Vec<String>,
    orientation: char,
    rotation_offset: usize,
    exact_position_matches: usize,
    candidates: Vec<crate::commands::resolve::RepeatPathConstraintCandidate>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PriorRepeatResolution {
    repeat_node: String,
    pairing_signature: String,
    source_candidate: String,
}

fn read_prior_repeat_resolution(
    paths: &OutputPaths,
    requested_repeat_node: &str,
) -> Option<PriorRepeatResolution> {
    let stats = fs::read_to_string(&paths.result_stats).ok()?;
    let mut repeat_node = None;
    let mut source_candidate = None;
    for line in stats.lines().skip(1) {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 5 || fields[0] != "summary" {
            continue;
        }
        match fields[3] {
            "repeat_resolution_node" => repeat_node = Some(fields[4].to_string()),
            "repeat_resolution_candidate" => source_candidate = Some(fields[4].to_string()),
            _ => {}
        }
    }
    let repeat_node = repeat_node?;
    let source_candidate = source_candidate?;
    if repeat_node != requested_repeat_node || source_candidate == "." {
        return None;
    }

    let resolution = fs::read_to_string(&paths.repeat_resolution).ok()?;
    let pairing_signature = resolution.lines().skip(1).find_map(|line| {
        let fields = line.split('\t').collect::<Vec<_>>();
        (fields.len() >= 3 && fields[0] == "summary" && fields[1] == repeat_node)
            .then(|| fields[2].to_string())
    })?;
    if pairing_signature == "." {
        return None;
    }

    Some(PriorRepeatResolution {
        repeat_node,
        pairing_signature,
        source_candidate,
    })
}

fn source_repeat_candidate(
    prior: Option<&PriorRepeatResolution>,
    decision: &RepeatResolutionDecision,
    local_candidate: &str,
) -> String {
    prior
        .filter(|prior| {
            prior.repeat_node == decision.repeat_node
                && prior.pairing_signature == decision.pairing_signature
        })
        .map(|prior| prior.source_candidate.clone())
        .unwrap_or_else(|| local_candidate.to_string())
}

pub(crate) fn select_dominant_repeat_pairing(
    rows: &[RepeatPathSupportRow],
    repeat_node: &str,
) -> Result<RepeatResolutionDecision, OrgraftError> {
    let target = rows
        .iter()
        .enumerate()
        .filter(|(_, row)| row.repeat_node == repeat_node)
        .collect::<Vec<_>>();
    if target.len() != 4 {
        return Err(OrgraftError::InvalidArgument(format!(
            "repeat `{repeat_node}` has {} complete path rows; expected 4 for a 2-in/2-out repeat",
            target.len()
        )));
    }

    let mut pairings = Vec::new();
    for left in 0..target.len() {
        for right in (left + 1)..target.len() {
            let (_, left_row) = target[left];
            let (_, right_row) = target[right];
            if left_row.left_edge_index == right_row.left_edge_index
                || left_row.right_edge_index == right_row.right_edge_index
            {
                continue;
            }
            let left_support = left_row.path_support.ok_or_else(|| {
                OrgraftError::InvalidArgument(format!(
                    "repeat `{repeat_node}` path support is unavailable; --reads is required"
                ))
            })?;
            let right_support = right_row.path_support.ok_or_else(|| {
                OrgraftError::InvalidArgument(format!(
                    "repeat `{repeat_node}` path support is unavailable; --reads is required"
                ))
            })?;
            let mut signatures = [
                repeat_path_physical_signature(left_row),
                repeat_path_physical_signature(right_row),
            ];
            signatures.sort();
            pairings.push((
                left_support + right_support,
                vec![format!("p{}", left + 1), format!("p{}", right + 1)],
                signatures.join("|"),
            ));
        }
    }
    pairings.sort_by(|left, right| right.0.cmp(&left.0).then_with(|| left.2.cmp(&right.2)));
    if pairings.len() != 2 {
        return Err(OrgraftError::InvalidArgument(format!(
            "repeat `{repeat_node}` produced {} perfect flank pairings; expected 2",
            pairings.len()
        )));
    }
    if pairings[0].0 == pairings[1].0 {
        return Err(OrgraftError::InvalidArgument(format!(
            "repeat `{repeat_node}` has tied complete-read pairing support ({})",
            pairings[0].0
        )));
    }

    Ok(RepeatResolutionDecision {
        repeat_node: repeat_node.to_string(),
        pairing_signature: pairings[0].2.clone(),
        selected_paths: pairings[0].1.clone(),
        selected_support: pairings[0].0,
        alternative_support: pairings[1].0,
    })
}

fn repeat_path_physical_signature(row: &RepeatPathSupportRow) -> String {
    let mut endpoints = [
        row.left_physical_endpoint.clone(),
        row.right_physical_endpoint.clone(),
    ];
    endpoints.sort();
    format!("{}--{}", endpoints[0], endpoints[1])
}

// pub(crate) + `debug_dir` instead of `&OutputPaths` so the sample-level
// refine module can reuse the identical four-path computation without
// a per-subgraph OutputPaths. Behaviour unchanged: keep_debug writes there,
// otherwise a temp work dir is used internally.
pub(crate) fn repeat_path_support_rows(
    gfa: &Gfa,
    node_classes: &BTreeMap<String, String>,
    reads: Option<&Path>,
    minimap2: &Path,
    threads: usize,
    debug_dir: &Path,
    keep_debug: bool,
) -> Result<Vec<RepeatPathSupportRow>, OrgraftError> {
    repeat_path_support_rows_with_extents(gfa, node_classes, &BTreeMap::new(), reads,
        minimap2, threads, debug_dir, keep_debug)
}

pub(crate) fn repeat_path_support_rows_with_extents(
    gfa: &Gfa,
    node_classes: &BTreeMap<String, String>,
    extents: &BTreeMap<String, (usize, usize)>,
    reads: Option<&Path>,
    minimap2: &Path,
    threads: usize,
    debug_dir: &Path,
    keep_debug: bool,
) -> Result<Vec<RepeatPathSupportRow>, OrgraftError> {
    repeat_path_support_rows_with_linear_chain_flanks(gfa, node_classes, extents, false,
        reads, minimap2, threads, debug_dir, keep_debug)
}

// Only refine's explicitly declared method enables chain flanks. The public
// rebuild entry and its existing callers retain adjacent-node templates.
pub(crate) fn repeat_path_support_rows_with_linear_chain_flanks(
    gfa: &Gfa, node_classes: &BTreeMap<String, String>,
    extents: &BTreeMap<String, (usize, usize)>, linear_chain_flanks: bool,
    reads: Option<&Path>, minimap2: &Path, threads: usize, debug_dir: &Path,
    keep_debug: bool,
) -> Result<Vec<RepeatPathSupportRow>, OrgraftError> {
    for node in extents.keys() {
        if !gfa.segments.contains_key(node) || node_classes.get(node).map(String::as_str) != Some("repeat_node") {
            return Err(OrgraftError::InvalidArgument("complete graph span names a missing or non-repeat node".into()));
        }
    }
    let mut rows = Vec::new();
    let mut templates = Vec::new();
    for (node, class) in node_classes {
        if class != "repeat_node" {
            continue;
        }
        let (left_extension_bp, right_extension_bp) = extents.get(node).copied().unwrap_or((0, 0));
        let full_span = left_extension_bp != 0 || right_extension_bp != 0;
        let mut left = Vec::new();
        let mut right = Vec::new();
        for (index, link) in gfa.links.iter().enumerate() {
            if let Some(side) = link_endpoint_side(link, node) {
                let item = repeat_path_edge_for_span(link, node, index + 1, full_span)?;
                if side == '-' {
                    left.push(item)
                } else {
                    right.push(item)
                }
            }
        }
        if full_span && (left.is_empty() || right.is_empty()) {
            return Err(OrgraftError::InvalidArgument("complete graph span lacks one of its actual outer arms".into()));
        }
        let status = if left.len() == 2 && right.len() == 2 {
            "unresolved"
        } else if left.len() <= 1 && right.len() <= 1 {
            "resolved"
        } else {
            "ambiguous"
        };
        if left.is_empty() || right.is_empty() {
            rows.push(RepeatPathSupportRow {
                repeat_node: node.clone(),
                repeat_status: status.to_string(),
                support_method: "unavailable_without_complete_path".to_string(),
                left_edge_index: None,
                right_edge_index: None,
                left_endpoint: ".".to_string(),
                right_endpoint: ".".to_string(),
                left_physical_endpoint: ".".to_string(),
                right_physical_endpoint: ".".to_string(),
                left_support: None,
                right_support: None,
                path_support: None,
                path_ratio: None,
                supporting_read_ids: Vec::new(),
            });
        }
        for left_edge in &left {
            for right_edge in &right {
                let row_index = rows.len();
                let id = format!("repeat_path_{:04}", templates.len() + 1);
                let (sequence, required_start, required_end) =
                    if linear_chain_flanks {
                        spell_repeat_path_linear_chain(gfa, node_classes, node, left_edge, right_edge,
                            left_extension_bp, right_extension_bp)?.0
                    } else {
                        spell_repeat_path_template_with_extents(gfa, node, left_edge, right_edge,
                            left_extension_bp, right_extension_bp)?
                    };
                templates.push(RepeatPathTemplate {
                    id,
                    repeat_node: node.clone(),
                    row_index,
                    sequence,
                    required_start,
                    required_end,
                });
                rows.push(RepeatPathSupportRow {
                    repeat_node: node.clone(),
                    repeat_status: status.to_string(),
                    support_method: if reads.is_some() {
                        if full_span { "spanning_reads_complete_graph_span".to_string() }
                        else if linear_chain_flanks { "spanning_reads_linear_chain_flanks".to_string() }
                        else { "spanning_reads".to_string() }
                    } else {
                        "unavailable_without_reads".to_string()
                    },
                    left_edge_index: Some(left_edge.index),
                    right_edge_index: Some(right_edge.index),
                    left_endpoint: left_edge.label(),
                    right_endpoint: right_edge.label(),
                    left_physical_endpoint: left_edge.physical_label(),
                    right_physical_endpoint: right_edge.physical_label(),
                    left_support: Some(left_edge.support),
                    right_support: Some(right_edge.support),
                    path_support: reads.map(|_| 0),
                    path_ratio: reads.map(|_| 0.0),
                    supporting_read_ids: Vec::new(),
                });
            }
        }
    }

    if let Some(reads) = reads {
        let support = map_reads_to_repeat_paths(
            &templates,
            reads,
            minimap2,
            threads,
            debug_dir,
            keep_debug,
        )?;
        for template in &templates {
            let mut read_ids = support.get(&template.id).cloned().unwrap_or_default();
            read_ids.sort();
            read_ids.dedup();
            rows[template.row_index].path_support = Some(read_ids.len());
            rows[template.row_index].supporting_read_ids = read_ids;
        }
        let mut totals = BTreeMap::<String, usize>::new();
        for row in &rows {
            if let Some(count) = row.path_support {
                *totals.entry(row.repeat_node.clone()).or_default() += count;
            }
        }
        for row in &mut rows {
            if let Some(count) = row.path_support {
                let total = totals.get(&row.repeat_node).copied().unwrap_or(0);
                row.path_ratio = Some(if total == 0 {
                    0.0
                } else {
                    count as f64 / total as f64
                });
            }
        }
    }
    Ok(rows)
}

fn repeat_path_edge(
    link: &Link,
    repeat_node: &str,
    index: usize,
) -> Result<RepeatPathEdge, OrgraftError> {
    let (endpoint_name, endpoint_orient) = if link.from_name == repeat_node {
        if link.from_orient == '+' {
            (link.to_name.clone(), link.to_orient)
        } else {
            (link.to_name.clone(), flip_orient(link.to_orient))
        }
    } else if link.to_name == repeat_node {
        if link.to_orient == '+' {
            (link.from_name.clone(), link.from_orient)
        } else {
            (link.from_name.clone(), flip_orient(link.from_orient))
        }
    } else {
        return Err(OrgraftError::InvalidArgument(format!(
            "link {index} does not touch repeat node `{repeat_node}`"
        )));
    };
    let endpoint_side = link_endpoint_side(link, &endpoint_name).ok_or_else(|| {
        OrgraftError::InvalidArgument(format!(
            "link {index} does not expose endpoint `{endpoint_name}`"
        ))
    })?;
    Ok(RepeatPathEdge {
        index,
        endpoint_name,
        endpoint_orient,
        endpoint_side,
        support: link_support(link),
        overlap: gfa_overlap_bases(&link.overlap),
    })
}

fn repeat_path_edge_for_span(link: &Link, node: &str, index: usize, full_span: bool)
    -> Result<RepeatPathEdge, OrgraftError>
{
    if !full_span { return repeat_path_edge(link, node, index); }
    let invalid = || OrgraftError::InvalidArgument("complete graph span has an unsupported or overflowing overlap".into());
    let mut length = 0usize;
    let mut number = 0usize;
    let mut digits = false;
    let mut operations = 0usize;
    for byte in link.overlap.bytes() {
        if byte.is_ascii_digit() {
            number = number.checked_mul(10).and_then(|value| value.checked_add((byte - b'0') as usize)).ok_or_else(invalid)?;
            digits = true;
        } else if matches!(byte, b'M' | b'=') && digits {
            length = length.checked_add(number).ok_or_else(invalid)?;
            number = 0;
            digits = false;
            operations = operations.checked_add(1).ok_or_else(invalid)?;
        } else { return Err(invalid()); }
    }
    if digits || operations == 0 { return Err(invalid()); }
    let mut orientation_only = link.clone();
    orientation_only.overlap = "0M".into();
    let mut edge = repeat_path_edge(&orientation_only, node, index)?;
    edge.overlap = length;
    Ok(edge)
}

fn spell_repeat_path_template_with_extents(
    gfa: &Gfa, repeat_node: &str, left: &RepeatPathEdge, right: &RepeatPathEdge,
    left_extension_bp: usize, right_extension_bp: usize,
) -> Result<(String, usize, usize), OrgraftError> {
    if left_extension_bp == 0 && right_extension_bp == 0 {
        return spell_repeat_path_template(gfa, repeat_node, left, right);
    }
    let invalid = || OrgraftError::InvalidArgument("complete graph span exceeds its exact arm or coordinate bounds".into());
    let left_sequence = oriented_segment_sequence(gfa, &left.endpoint_name, left.endpoint_orient)?;
    let repeat_sequence = oriented_segment_sequence(gfa, repeat_node, '+')?;
    let right_sequence = oriented_segment_sequence(gfa, &right.endpoint_name, right.endpoint_orient)?;
    if !left_sequence.is_ascii() || !repeat_sequence.is_ascii() || !right_sequence.is_ascii()
        || repeat_sequence.is_empty() || left.overlap > left_sequence.len() || left.overlap > repeat_sequence.len()
        || right.overlap > right_sequence.len() || right.overlap > repeat_sequence.len()
    { return Err(invalid()); }
    if (left.overlap != 0 && left_sequence[left_sequence.len()-left.overlap..] != repeat_sequence[..left.overlap])
        || (right.overlap != 0 && repeat_sequence[repeat_sequence.len()-right.overlap..] != right_sequence[..right.overlap])
    { return Err(OrgraftError::InvalidArgument("complete graph span has inconsistent actual overlap DNA".into())); }
    let left_available = left_sequence.len().checked_sub(left.overlap).ok_or_else(invalid)?;
    let right_available = right_sequence.len().checked_sub(right.overlap).ok_or_else(invalid)?;
    let left_needed = left_extension_bp.checked_add(REPEAT_PATH_REQUIRED_FLANK_BP).ok_or_else(invalid)?;
    let right_needed = right_extension_bp.checked_add(REPEAT_PATH_REQUIRED_FLANK_BP).ok_or_else(invalid)?;
    if left_available < left_needed || right_available < right_needed {
        return Err(OrgraftError::InvalidArgument("complete graph span lacks an extension plus the full 1000-bp outer anchor".into()));
    }
    let left_take = left_extension_bp.checked_add(REPEAT_PATH_TEMPLATE_FLANK_BP)
        .and_then(|value| value.checked_add(left.overlap)).ok_or_else(invalid)?.min(left_sequence.len());
    let right_take = right_extension_bp.checked_add(REPEAT_PATH_TEMPLATE_FLANK_BP)
        .and_then(|value| value.checked_add(right.overlap)).ok_or_else(invalid)?.min(right_sequence.len());
    let left_start = left_sequence.len().checked_sub(left_take).ok_or_else(invalid)?;
    let repeat_start = left_take.checked_sub(left.overlap).ok_or_else(invalid)?;
    let repeat_end = repeat_start.checked_add(repeat_sequence.len()).ok_or_else(invalid)?;
    let full_start = repeat_start.checked_sub(left_extension_bp).ok_or_else(invalid)?;
    let full_end = repeat_end.checked_add(right_extension_bp).ok_or_else(invalid)?;
    let required_start = full_start.checked_sub(REPEAT_PATH_REQUIRED_FLANK_BP).ok_or_else(invalid)?;
    let required_end = full_end.checked_add(REPEAT_PATH_REQUIRED_FLANK_BP).ok_or_else(invalid)?;
    let capacity = left_take.checked_add(repeat_sequence.len().checked_sub(left.overlap).ok_or_else(invalid)?)
        .and_then(|value| value.checked_add(right_take.checked_sub(right.overlap)?)).ok_or_else(invalid)?;
    if required_end > capacity { return Err(invalid()); }
    let mut sequence = String::with_capacity(capacity);
    sequence.push_str(&left_sequence[left_start..]);
    sequence.push_str(&repeat_sequence[left.overlap..]);
    sequence.push_str(&right_sequence[right.overlap..right_take]);
    Ok((sequence, required_start, required_end))
}

fn spell_repeat_path_template(
    gfa: &Gfa,
    repeat_node: &str,
    left: &RepeatPathEdge,
    right: &RepeatPathEdge,
) -> Result<(String, usize, usize), OrgraftError> {
    let left_sequence = oriented_segment_sequence(gfa, &left.endpoint_name, left.endpoint_orient)?;
    let repeat_sequence = oriented_segment_sequence(gfa, repeat_node, '+')?;
    let right_sequence =
        oriented_segment_sequence(gfa, &right.endpoint_name, right.endpoint_orient)?;
    let left_overlap = left
        .overlap
        .min(left_sequence.len())
        .min(repeat_sequence.len());
    let right_overlap = right
        .overlap
        .min(right_sequence.len())
        .min(repeat_sequence.len());
    let left_take = (REPEAT_PATH_TEMPLATE_FLANK_BP + left_overlap).min(left_sequence.len());
    let right_take = (REPEAT_PATH_TEMPLATE_FLANK_BP + right_overlap).min(right_sequence.len());
    let left_chunk = &left_sequence[left_sequence.len() - left_take..];
    let right_chunk = &right_sequence[..right_take];

    let mut sequence =
        String::with_capacity(left_chunk.len() + repeat_sequence.len() + right_chunk.len());
    sequence.push_str(left_chunk);
    sequence.push_str(&repeat_sequence[left_overlap..]);
    sequence.push_str(&right_chunk[right_overlap..]);

    let repeat_start = left_chunk.len().saturating_sub(left_overlap);
    let repeat_end = left_chunk.len() + repeat_sequence.len().saturating_sub(left_overlap);
    let left_required = REPEAT_PATH_REQUIRED_FLANK_BP.min(repeat_start);
    let right_required = REPEAT_PATH_REQUIRED_FLANK_BP.min(sequence.len() - repeat_end);
    Ok((
        sequence,
        repeat_start - left_required,
        repeat_end + right_required,
    ))
}

pub(crate) fn repeat_node_path_templates(gfa: &Gfa, node: &str)
    -> Result<Vec<(usize, usize, String, usize, usize)>, OrgraftError>
{
    repeat_node_path_templates_with_extents(gfa, node, 0, 0)
}

pub(crate) fn repeat_node_path_templates_with_extents(
    gfa: &Gfa, node: &str, left_extension_bp: usize, right_extension_bp: usize,
) -> Result<Vec<(usize, usize, String, usize, usize)>, OrgraftError> {
    let full_span = left_extension_bp != 0 || right_extension_bp != 0;
    let mut left = Vec::new();
    let mut right = Vec::new();
    for (index, link) in gfa.links.iter().enumerate() {
        if let Some(side) = link_endpoint_side(link, node) {
            let edge = repeat_path_edge_for_span(link, node, index + 1, full_span)?;
            if side == '-' { left.push(edge); } else { right.push(edge); }
        }
    }
    if left.len() != 2 || right.len() != 2 {
        return Err(OrgraftError::InvalidArgument("template proof requires four physical 2x2 paths".into()));
    }
    let mut result = Vec::new();
    for a in &left { for b in &right {
        let (sequence, start, end) = spell_repeat_path_template_with_extents(gfa, node, a, b, left_extension_bp, right_extension_bp)?;
        result.push((a.index, b.index, sequence, start, end));
    }}
    Ok(result)
}


// Follow an outward arm only as far as the existing template length requires.
// Every traversed node must be a single occurrence in the verified source walk.
fn linear_chain_arm(
    gfa: &Gfa, classes: &BTreeMap<String, String>, name: &str, orient: char,
    needed: usize,
) -> Result<(String, String), OrgraftError> {
    let invalid = || OrgraftError::InvalidArgument(
        "linear-chain flank lacks a unique source-bound non-repeat 0M continuation".into());
    let mut current = (name.to_string(), orient);
    let mut seen = BTreeSet::new();
    let mut sequence = String::new();
    let mut path = Vec::new();
    loop {
        if !seen.insert(current.0.clone())
            || classes.get(&current.0).map(String::as_str) != Some("single_copy_node")
            || !matches!(current.1, '+' | '-') { return Err(invalid()); }
        let dna = oriented_segment_sequence(gfa, &current.0, current.1)?;
        if dna.is_empty() || !dna.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A'|b'C'|b'G'|b'T')) {
            return Err(invalid());
        }
        sequence.push_str(&dna);
        path.push(format!("{}{}", current.0, current.1));
        if sequence.len() >= needed { break; }
        let outgoing = |node: &(String, char)| {
            let mut edges = Vec::new();
            for (i, link) in gfa.links.iter().enumerate() {
                if (link.from_name.as_str(), link.from_orient) == (node.0.as_str(), node.1) {
                    edges.push((i, (link.to_name.clone(), link.to_orient)));
                }
                if (link.to_name.as_str(), flip_orient(link.to_orient)) == (node.0.as_str(), node.1) {
                    edges.push((i, (link.from_name.clone(), flip_orient(link.from_orient))));
                }
            }
            edges
        };
        if outgoing(&(current.0.clone(), flip_orient(current.1))).len() != 1 {
            return Err(invalid());
        }
        let next = outgoing(&current);
        if next.len() != 1 || gfa.links[next[0].0].overlap != "0M" { return Err(invalid()); }
        let incoming = outgoing(&(next[0].1.0.clone(), flip_orient(next[0].1.1)));
        if incoming.len() != 1 || incoming[0].0 != next[0].0 { return Err(invalid()); }
        path.push(format!("edge{}", next[0].0 + 1));
        current = next[0].1.clone();
    }
    Ok((sequence, path.join(",")))
}

fn spell_repeat_path_linear_chain(
    gfa: &Gfa, classes: &BTreeMap<String, String>, node: &str,
    left: &RepeatPathEdge, right: &RepeatPathEdge, le: usize, re: usize,
) -> Result<((String, usize, usize), String), OrgraftError> {
    let invalid = || OrgraftError::InvalidArgument("linear-chain template has invalid overlap or bounds".into());
    // New traversal is confined to zero-overlap arms; it never normalizes a graph.
    if gfa.links[left.index-1].overlap != "0M" || gfa.links[right.index-1].overlap != "0M" {
        return Err(invalid());
    }
    let ln = le.checked_add(REPEAT_PATH_TEMPLATE_FLANK_BP).ok_or_else(invalid)?;
    let rn = re.checked_add(REPEAT_PATH_TEMPLATE_FLANK_BP).ok_or_else(invalid)?;
    let (left_out, lp) = linear_chain_arm(gfa, classes, &left.endpoint_name, flip_orient(left.endpoint_orient), ln)?;
    let (right_out, rp) = linear_chain_arm(gfa, classes, &right.endpoint_name, right.endpoint_orient, rn)?;
    let left_word = reverse_complement(&left_out[..ln]);
    let repeat = oriented_segment_sequence(gfa, node, '+')?;
    if repeat.is_empty() { return Err(invalid()); }
    let end = ln.checked_add(repeat.len()).and_then(|n| n.checked_add(re)).and_then(|n| n.checked_add(REPEAT_PATH_REQUIRED_FLANK_BP)).ok_or_else(invalid)?;
    let word = format!("{left_word}{repeat}{}", &right_out[..rn]);
    Ok(((word, REPEAT_PATH_TEMPLATE_FLANK_BP-REPEAT_PATH_REQUIRED_FLANK_BP, end),
        format!("left={lp};right={rp}")))
}

pub(crate) fn repeat_node_linear_chain_templates(
    gfa: &Gfa, classes: &BTreeMap<String, String>, node: &str, le: usize, re: usize,
) -> Result<Vec<(usize, usize, String, usize, usize, String)>, OrgraftError> {
    let mut left = Vec::new();
    let mut right = Vec::new();
    for (i, link) in gfa.links.iter().enumerate() {
        if let Some(side) = link_endpoint_side(link, node) {
            let edge = repeat_path_edge_for_span(link, node, i+1, true)?;
            if side == '-' { left.push(edge); } else { right.push(edge); }
        }
    }
    if classes.get(node).map(String::as_str) != Some("repeat_node") || left.len()!=2 || right.len()!=2 {
        return Err(OrgraftError::InvalidArgument("linear-chain proof requires a source-bound repeat and all four physical paths".into()));
    }
    let mut paths = Vec::new();
    for a in &left { for b in &right {
        let ((dna,start,end),proof) = spell_repeat_path_linear_chain(gfa, classes, node, a, b, le, re)?;
        paths.push((a.index,b.index,dna,start,end,proof));
    }}
    Ok(paths)
}

fn oriented_segment_sequence(gfa: &Gfa, name: &str, orient: char) -> Result<String, OrgraftError> {
    let segment = gfa.segments.get(name).ok_or_else(|| {
        OrgraftError::InvalidArgument(format!("repeat path references missing segment `{name}`"))
    })?;
    Ok(if orient == '+' {
        segment.sequence.clone()
    } else {
        reverse_complement(&segment.sequence)
    })
}

fn gfa_overlap_bases(value: &str) -> usize {
    let mut number = String::new();
    let mut overlap = 0usize;
    for character in value.chars() {
        if character.is_ascii_digit() {
            number.push(character);
            continue;
        }
        let length = number.parse::<usize>().unwrap_or(0);
        if matches!(character, 'M' | '=' | 'X') {
            overlap += length;
        }
        number.clear();
    }
    overlap
}

#[derive(Clone, Debug)]
pub(crate) struct RepeatPathContext {
    pub(crate) left_edge_index: usize,
    pub(crate) right_edge_index: usize,
    pub(crate) sequence: String,
    pub(crate) required_start: usize,
    pub(crate) required_end: usize,
    pub(crate) source_path: String,
}

// Enumerate every physical continuation until the unchanged native template
// flank is available. A limit or any incomplete branch rejects the whole node.
pub(crate) fn repeat_node_path_contexts(
    gfa: &Gfa,
    node: &str,
    occurrences: &BTreeMap<String, usize>,
    left_extension: usize,
    right_extension: usize,
) -> Result<Vec<RepeatPathContext>, OrgraftError> {
    let invalid = || {
        OrgraftError::InvalidArgument(
            "multi-context flank lacks complete source-bound physical continuations".into(),
        )
    };
    fn outgoing(gfa: &Gfa, n: &str, s: char) -> Vec<(usize, String, char)> {
        let mut v = Vec::new();
        for (i, l) in gfa.links.iter().enumerate() {
            if l.from_name == n && l.from_orient == s {
                v.push((i, l.to_name.clone(), l.to_orient));
            }
            if l.to_name == n && flip_orient(l.to_orient) == s {
                v.push((i, l.from_name.clone(), flip_orient(l.from_orient)));
            }
        }
        v
    }
    fn arm(
        gfa: &Gfa,
        occurrences: &BTreeMap<String, usize>,
        item: (usize, String, char),
        word: String,
        path: Vec<String>,
        edges: BTreeSet<usize>,
        seen: BTreeSet<String>,
        needed: usize,
        result: &mut Vec<(String, String)>,
    ) -> Result<(), OrgraftError> {
        let invalid = || {
            OrgraftError::InvalidArgument(
                "multi-context path is cyclic, incomplete, unbound, or exceeds64 contexts".into(),
            )
        };
        let (i, n, s) = item;
        if edges.contains(&i)
            || seen.contains(&n)
            || gfa.links[i].overlap != "0M"
            || seen.len() > 64
        {
            return Err(invalid());
        }
        let degree = (outgoing(gfa, &n, '-').len(), outgoing(gfa, &n, '+').len());
        if !matches!(
            (occurrences.get(&n).copied(), degree),
            (Some(1), (1, 1)) | (Some(2), (2, 2))
        ) {
            return Err(invalid());
        }
        let dna = oriented_segment_sequence(gfa, &n, s)?;
        if dna.is_empty() || !dna.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            return Err(invalid());
        }
        let mut word = word;
        word.push_str(&dna);
        let mut path = path;
        path.push(format!("{}:{}{}", i + 1, n, s));
        if word.len() >= needed {
            if result.len() >= 64 {
                return Err(invalid());
            }
            result.push((word[..needed].to_string(), path.join(",")));
            return Ok(());
        }
        let mut edges = edges;
        edges.insert(i);
        let mut seen = seen;
        seen.insert(n.clone());
        let next = outgoing(gfa, &n, s);
        if next.is_empty() {
            return Err(invalid());
        }
        for item in next {
            arm(
                gfa,
                occurrences,
                item,
                word.clone(),
                path.clone(),
                edges.clone(),
                seen.clone(),
                needed,
                result,
            )?;
        }
        Ok(())
    }
    if occurrences.get(node) != Some(&2) {
        return Err(invalid());
    }
    let core = oriented_segment_sequence(gfa, node, '+')?;
    if core.is_empty() || !core.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
        return Err(invalid());
    }
    let left = outgoing(gfa, node, '-');
    let right = outgoing(gfa, node, '+');
    if left.len() != 2 || right.len() != 2 {
        return Err(invalid());
    }
    let mut result = Vec::new();
    let mut unique = BTreeSet::new();
    for a in left {
        for b in &right {
            let mut ls = Vec::new();
            let mut rs = Vec::new();
            arm(
                gfa,
                occurrences,
                a.clone(),
                String::new(),
                Vec::new(),
                BTreeSet::new(),
                BTreeSet::from([node.to_string()]),
                REPEAT_PATH_TEMPLATE_FLANK_BP + left_extension,
                &mut ls,
            )?;
            arm(
                gfa,
                occurrences,
                b.clone(),
                String::new(),
                Vec::new(),
                BTreeSet::new(),
                BTreeSet::from([node.to_string()]),
                REPEAT_PATH_TEMPLATE_FLANK_BP + right_extension,
                &mut rs,
            )?;
            for (l, lp) in &ls {
                for (r, rp) in &rs {
                    if result.len() >= 256 {
                        return Err(invalid());
                    }
                    let sequence = format!("{}{}{}", reverse_complement(l), core, r);
                    // Distinct physical contexts must be distinguishable in full DNA.
                    if !unique.insert(sequence.clone()) {
                        return Err(invalid());
                    }
                    result.push(RepeatPathContext {
                        left_edge_index: a.0 + 1,
                        right_edge_index: b.0 + 1,
                        required_start: REPEAT_PATH_TEMPLATE_FLANK_BP
                            - REPEAT_PATH_REQUIRED_FLANK_BP,
                        required_end: REPEAT_PATH_TEMPLATE_FLANK_BP
                            + left_extension
                            + core.len()
                            + right_extension
                            + REPEAT_PATH_REQUIRED_FLANK_BP,
                        sequence,
                        source_path: format!("left:{lp}|right:{rp}"),
                    });
                }
            }
        }
    }
    Ok(result)
}

pub(crate) fn repeat_path_support_rows_with_contexts(
    gfa: &Gfa,
    classes: &BTreeMap<String, String>,
    occurrences: &BTreeMap<String, usize>,
    extents: &BTreeMap<String, (usize, usize)>,
    reads: &Path,
    minimap2: &Path,
    threads: usize,
    debug_dir: &Path,
    keep_debug: bool,
) -> Result<Vec<RepeatPathSupportRow>, OrgraftError> {
    let mut rows = repeat_path_support_rows_with_extents(
        gfa,
        classes,
        &BTreeMap::new(),
        None,
        minimap2,
        threads,
        debug_dir,
        keep_debug,
    )?;
    let mut templates = Vec::new();
    for (node, class) in classes {
        if class == "repeat_node" {
            let (le, re) = extents.get(node).copied().unwrap_or((0, 0));
            for c in repeat_node_path_contexts(gfa, node, occurrences, le, re)? {
                let indices = rows
                    .iter()
                    .enumerate()
                    .filter(|(_, r)| {
                        r.repeat_node == *node
                            && r.left_edge_index == Some(c.left_edge_index)
                            && r.right_edge_index == Some(c.right_edge_index)
                    })
                    .map(|(i, _)| i)
                    .collect::<Vec<_>>();
                if indices.len() != 1 {
                    return Err(OrgraftError::InvalidArgument(
                        "context lacks one physical port row".into(),
                    ));
                }
                templates.push(RepeatPathTemplate {
                    id: format!("repeat_path_{:04}", templates.len() + 1),
                    repeat_node: node.clone(),
                    row_index: indices[0],
                    sequence: c.sequence,
                    required_start: c.required_start,
                    required_end: c.required_end,
                });
            }
        }
    }
    let support = map_reads_to_repeat_paths_grouped(
        &templates, reads, minimap2, threads, debug_dir, keep_debug, true,
    )?;
    let mut union = BTreeMap::<usize, BTreeSet<String>>::new();
    for t in &templates {
        union
            .entry(t.row_index)
            .or_default()
            .extend(support.get(&t.id).cloned().unwrap_or_default());
    }
    let mut totals = BTreeMap::<String, usize>::new();
    for (i, r) in rows.iter_mut().enumerate() {
        r.supporting_read_ids = union.remove(&i).unwrap_or_default().into_iter().collect();
        r.path_support = Some(r.supporting_read_ids.len());
        r.support_method = "spanning_reads_complete_physical_contexts".into();
        *totals.entry(r.repeat_node.clone()).or_default() += r.supporting_read_ids.len();
    }
    for r in &mut rows {
        let n = totals[&r.repeat_node];
        r.path_ratio = Some(if n == 0 {
            0.0
        } else {
            r.supporting_read_ids.len() as f64 / n as f64
        });
    }
    Ok(rows)
}

fn map_reads_to_repeat_paths(
    templates: &[RepeatPathTemplate],
    reads: &Path,
    minimap2: &Path,
    threads: usize,
    debug_dir: &Path,
    keep_debug: bool,
) -> Result<BTreeMap<String, Vec<String>>, OrgraftError> {
    map_reads_to_repeat_paths_grouped(templates,reads,minimap2,threads,debug_dir,keep_debug,false)
}

fn map_reads_to_repeat_paths_grouped(
    templates: &[RepeatPathTemplate],
    reads: &Path,
    minimap2: &Path,
    threads: usize,
    debug_dir: &Path,
    keep_debug: bool,
    physical_contexts: bool,
) -> Result<BTreeMap<String, Vec<String>>, OrgraftError> {
    if templates.is_empty() {
        return Ok(BTreeMap::new());
    }
    let temp_dir = if keep_debug {
        debug_dir.to_path_buf()
    } else {
        temp_work_dir("orgraft-repeat-path-support")
    };
    fs::create_dir_all(&temp_dir)?;
    let targets = temp_dir.join("repeat_path_targets.fasta");
    let mut target_file = File::create(&targets)?;
    for template in templates {
        write_fasta_record(&mut target_file, &template.id, &template.sequence)?;
    }
    let output = Command::new(minimap2)
        .args([
            "-x",
            "map-hifi",
            "-c",
            "--eqx",
            "--secondary=yes",
            "-N",
            &REPEAT_PATH_SECONDARY_LIMIT.to_string(),
            "-t",
            &threads.max(1).to_string(),
        ])
        .arg(&targets)
        .arg(reads)
        .output()
        .map_err(|error| {
            OrgraftError::InvalidArgument(format!(
                "failed to map reads for repeat path support: {error}"
            ))
        })?;
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();
    if keep_debug {
        fs::write(debug_dir.join("repeat_path_support.paf"), &stdout)?;
        fs::write(
            debug_dir.join("repeat_path_support.stderr.log"),
            &stderr,
        )?;
    }
    if !output.status.success() {
        if !keep_debug {
            let _ = fs::remove_dir_all(&temp_dir);
        }
        return Err(OrgraftError::InvalidArgument(format!(
            "minimap2 repeat path support mapping failed: {stderr}"
        )));
    }

    let support = repeat_path_support_from_paf(templates, &stdout, physical_contexts);
    if !keep_debug { let _ = fs::remove_dir_all(&temp_dir); }
    support
}

// The existing alignment filter, ranking and ambiguity margins are shared.
// Contexts for one physical port pair compete once, never as extra paths.
fn repeat_path_support_from_paf(templates: &[RepeatPathTemplate], stdout: &str,
    physical_contexts: bool) -> Result<BTreeMap<String,Vec<String>>,OrgraftError> {
    // A saturated secondary list cannot prove that every competing physical
    // row was emitted. Reject the entire context measurement; do not discard
    // that read and turn a censored query into a lower support count.
    if physical_contexts {
        let mut secondary_counts = BTreeMap::<String, usize>::new();
        for line in stdout.lines().filter(|line| !line.trim().is_empty()) {
            let alignment = parse_repeat_path_paf(line)?;
            if !alignment.primary {
                let count = secondary_counts.entry(alignment.read_id.clone()).or_default();
                *count += 1;
                if *count >= REPEAT_PATH_SECONDARY_LIMIT {
                    return Err(OrgraftError::InvalidArgument(format!(
                        "context measurement has a saturated secondary alignment list for read {}; physical-port competition is not complete",
                        alignment.read_id
                    )));
                }
            }
        }
    }
    let by_id = templates
        .iter()
        .map(|template| (template.id.as_str(), template))
        .collect::<HashMap<_, _>>();
    let mut by_read_repeat = BTreeMap::<(String, String), Vec<RepeatPathAlignment>>::new();
    for line in stdout.lines().filter(|line| !line.trim().is_empty()) {
        let alignment = parse_repeat_path_paf(line)?;
        let Some(template) = by_id.get(alignment.target_id.as_str()) else {
            continue;
        };
        let tolerance = REPEAT_PATH_ALIGNMENT_TOLERANCE_BP;
        if alignment.identity() < REPEAT_PATH_MIN_IDENTITY
            || alignment.target_start > template.required_start.saturating_add(tolerance)
            || alignment.target_end.saturating_add(tolerance) < template.required_end
        {
            continue;
        }
        by_read_repeat
            .entry((alignment.read_id.clone(), template.repeat_node.clone()))
            .or_default()
            .push(alignment);
    }

    let mut support = BTreeMap::<String, Vec<String>>::new();
    for ((read_id, _), alignments) in by_read_repeat {
        let mut best_by_target = BTreeMap::<String, RepeatPathAlignment>::new();
        for alignment in alignments {
            let key = if physical_contexts {format!("row:{}",by_id[alignment.target_id.as_str()].row_index)} else {alignment.target_id.clone()};
            let replace = best_by_target
                .get(&key)
                .is_none_or(|old| repeat_path_alignment_better(&alignment, old));
            if replace {
                best_by_target.insert(key, alignment);
            }
        }
        let mut candidates = best_by_target.into_values().collect::<Vec<_>>();
        candidates.sort_by(|left, right| {
            repeat_path_alignment_order(right, left)
                .then_with(|| left.target_id.cmp(&right.target_id))
        });
        let Some(best) = candidates.first() else {
            continue;
        };
        if candidates.get(1).is_some_and(|second| {
            best.matches.abs_diff(second.matches) <= 20
                && best.target_span().abs_diff(second.target_span()) <= 50
        }) {
            continue;
        }
        support
            .entry(best.target_id.clone())
            .or_default()
            .push(read_id);
    }
    Ok(support)
}

fn parse_repeat_path_paf(line: &str) -> Result<RepeatPathAlignment, OrgraftError> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() < 12 {
        return Err(OrgraftError::InvalidArgument(format!(
            "malformed repeat path PAF line: {line}"
        )));
    }
    Ok(RepeatPathAlignment {
        read_id: fields[0].to_string(),
        target_id: fields[5].to_string(),
        query_start: parse_usize_value(fields[2], "PAF query start")?,
        query_end: parse_usize_value(fields[3], "PAF query end")?,
        target_start: parse_usize_value(fields[7], "PAF target start")?,
        target_end: parse_usize_value(fields[8], "PAF target end")?,
        matches: parse_usize_value(fields[9], "PAF matches")?,
        mapq: parse_usize_value(fields[11], "PAF mapq")?,
        primary: fields.iter().any(|field| *field == "tp:A:P"),
    })
}

fn repeat_path_alignment_better(left: &RepeatPathAlignment, right: &RepeatPathAlignment) -> bool {
    repeat_path_alignment_order(left, right) == Ordering::Greater
}

fn repeat_path_alignment_order(
    left: &RepeatPathAlignment,
    right: &RepeatPathAlignment,
) -> Ordering {
    left.matches
        .cmp(&right.matches)
        .then_with(|| left.target_span().cmp(&right.target_span()))
        .then_with(|| left.primary.cmp(&right.primary))
        .then_with(|| left.mapq.cmp(&right.mapq))
}

fn write_repeat_path_support(
    path: &Path,
    rows: &[RepeatPathSupportRow],
) -> Result<(), OrgraftError> {
    let mut output = String::from(
        "repeat_node\tpath\tleft_endpoint\tright_endpoint\tmethod\tspanning_reads\tpath_ratio\tread_ids\n",
    );
    let mut counters = BTreeMap::<&str, usize>::new();
    for row in rows {
        let counter = counters.entry(&row.repeat_node).or_default();
        *counter += 1;
        writeln!(
            output,
            "{}\tp{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.repeat_node,
            counter,
            row.left_endpoint,
            row.right_endpoint,
            row.support_method,
            optional_usize(row.path_support),
            optional_f64(row.path_ratio),
            if row.supporting_read_ids.is_empty() {
                ".".to_string()
            } else {
                row.supporting_read_ids.join(",")
            },
        )
        .unwrap();
    }
    fs::write(path, output)?;
    Ok(())
}

fn write_repeat_resolution(
    path: &Path,
    outcome: Option<&RepeatResolutionOutcome>,
) -> Result<(), OrgraftError> {
    let mut output = String::from(
        "section\titem\tpairing_signature\tselected\tscore\tcontinuous_fraction\tdetails\n",
    );
    let Some(outcome) = outcome else {
        output.push_str("summary\tdisabled\t.\tNo\t.\t.\trepeat resolution not requested\n");
        fs::write(path, output)?;
        return Ok(());
    };
    writeln!(
        output,
        "summary\t{}\t{}\tYes\t{}\t{}\tmethod={};source_candidate={};local_candidate={};selected_paths={};selected_support={};alternative_support={};orientation={};rotation_offset={};exact_position_matches={};selected_order={};ordered_nodes={};oriented_nodes={}",
        outcome.decision.repeat_node,
        outcome.decision.pairing_signature,
        optional_f64(outcome.score),
        optional_f64(outcome.continuous_fraction),
        outcome.selection_method,
        outcome.selected_candidate,
        outcome.local_candidate,
        outcome.decision.selected_paths.join(","),
        outcome.decision.selected_support,
        outcome.decision.alternative_support,
        outcome.orientation,
        outcome.rotation_offset,
        outcome.exact_position_matches,
        outcome.selected_order,
        outcome.ordered_nodes.join(","),
        outcome.oriented_nodes.join(","),
    )
    .unwrap();
    for candidate in &outcome.candidates {
        writeln!(
            output,
            "candidate\t{}\t{}\t{}\t{}\t{}\tfiltered_candidate_count={};total_candidate_count={}",
            candidate.candidate,
            candidate.pairing_signature,
            if candidate.selected { "Yes" } else { "No" },
            optional_f64(candidate.score),
            optional_f64(candidate.continuous_fraction),
            outcome.filtered_candidate_count,
            outcome.candidate_count,
        )
        .unwrap();
    }
    fs::write(path, output)?;
    Ok(())
}

fn attach_repeat_path_support_paths(gfa: &mut Gfa, rows: &[RepeatPathSupportRow]) {
    gfa.other_lines
        .retain(|fields| fields.first().map(String::as_str) != Some("P"));
    let mut methods = BTreeMap::<&str, &str>::new();
    for row in rows {
        methods.insert(
            &row.repeat_node,
            if row.path_support.is_some() {
                "spanning_reads"
            } else {
                "unavailable"
            },
        );
    }
    for (node, method) in methods {
        if let Some(segment) = gfa.segments.get_mut(node) {
            segment.tags = replace_tag(&segment.tags, format!("PM:Z:{method}"));
        }
    }

    let mut counters: BTreeMap<&str, usize> = BTreeMap::new();
    for row in rows {
        let (Some(_), Some(_), Some(path_support)) =
            (row.left_edge_index, row.right_edge_index, row.path_support)
        else {
            continue;
        };
        let counter = counters.entry(&row.repeat_node).or_default();
        *counter += 1;
        let path_index = format!("p{counter}");
        let path = format!(
            "{},{},{}",
            row.left_endpoint,
            oriented_node(&row.repeat_node, '+'),
            row.right_endpoint
        );
        gfa.other_lines.push(vec![
            "P".to_string(),
            format!("repeat_{}_{}", row.repeat_node, path_index),
            path,
            "*,*".to_string(),
            "PT:Z:repeat_path_support".to_string(),
            format!("RN:Z:{}", row.repeat_node),
            format!("PI:Z:{path_index}"),
            format!("RS:Z:{}", row.repeat_status),
            "PM:Z:spanning_reads".to_string(),
            format!("RC:i:{path_support}"),
            format!("PR:f:{:.6}", row.path_ratio.unwrap_or(0.0)),
            format!("LE:Z:{}", row.left_endpoint),
            format!("RE:Z:{}", row.right_endpoint),
        ]);
    }
}

pub(crate) fn oriented_node(name: &str, orient: char) -> String {
    format!("{name}{orient}")
}

fn optional_usize(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| ".".to_string())
}

fn optional_f64(value: Option<f64>) -> String {
    value
        .map(|value| format!("{value:.6}"))
        .unwrap_or_else(|| ".".to_string())
}

#[derive(Debug)]
struct CoordinateConsistency {
    polished_length: usize,
    polished_hash64: String,
    mapped_segments: usize,
    unmapped_segments: usize,
    covered_bases: usize,
    multi_covered_bases: usize,
    gap_bases: usize,
    coverage_fraction: f64,
    raw_segments: usize,
    merged_segments: usize,
    verified_segments: usize,
}

impl CoordinateConsistency {
    fn new(
        record: &FastaRecord,
        raw: &Gfa,
        merged: &Gfa,
        verified: &Gfa,
        mapping: &Mapping,
    ) -> Self {
        let mut coverage = vec![0usize; record.sequence.len()];
        let mut mapped = 0usize;
        let mut unmapped = 0usize;
        for name in &merged.order {
            let node_mapping = &mapping.by_node[name];
            let hits = if node_mapping.accepted_hits.is_empty() {
                node_mapping.selected.iter().collect::<Vec<_>>()
            } else {
                node_mapping.accepted_hits.iter().collect::<Vec<_>>()
            };
            if hits.is_empty() {
                unmapped += 1;
            } else {
                mapped += 1;
            }
            for hit in hits {
                for pos in hit.target_start..hit.target_end {
                    coverage[pos % record.sequence.len()] += 1;
                }
            }
        }
        let covered = coverage.iter().filter(|value| **value > 0).count();
        let multi = coverage.iter().filter(|value| **value > 1).count();
        let gap_bases = record.sequence.len().saturating_sub(covered);
        Self {
            polished_length: record.sequence.len(),
            polished_hash64: stable_hash64(&record.sequence),
            mapped_segments: mapped,
            unmapped_segments: unmapped,
            covered_bases: covered,
            multi_covered_bases: multi,
            gap_bases,
            coverage_fraction: covered as f64 / record.sequence.len() as f64,
            raw_segments: raw.order.len(),
            merged_segments: merged.order.len(),
            verified_segments: verified.order.len(),
        }
    }

    fn linear_tiling_status(&self) -> &'static str {
        if self.gap_bases == 0 && self.multi_covered_bases == 0 {
            "PASS"
        } else {
            "WARN"
        }
    }
}

#[derive(Debug)]
struct GraphStats {
    segments: usize,
    links: usize,
    bp: usize,
}

fn graph_stats(gfa: &Gfa) -> GraphStats {
    GraphStats {
        segments: gfa.order.len(),
        links: gfa.links.len(),
        bp: gfa
            .order
            .iter()
            .map(|name| gfa.segments[name].sequence.len())
            .sum(),
    }
}

fn parse_eqx_counts(cigar: &str) -> Result<(usize, usize), OrgraftError> {
    let mut length = 0usize;
    let mut has_length = false;
    let mut exact = 0usize;
    let mut mismatch = 0usize;
    for character in cigar.chars() {
        if character.is_ascii_digit() {
            length = length
                .checked_mul(10)
                .and_then(|value| value.checked_add(character.to_digit(10).unwrap() as usize))
                .ok_or_else(|| {
                    OrgraftError::InvalidArgument(format!("CIGAR length overflow: {cigar}"))
                })?;
            has_length = true;
            continue;
        }
        if !has_length {
            return Err(OrgraftError::InvalidArgument(format!(
                "malformed CIGAR operation in `{cigar}`"
            )));
        }
        match character {
            '=' => exact = exact.saturating_add(length),
            'X' => mismatch = mismatch.saturating_add(length),
            'M' | 'I' | 'D' | 'N' | 'S' | 'H' | 'P' => {}
            _ => {
                return Err(OrgraftError::InvalidArgument(format!(
                    "unsupported CIGAR operation `{character}` in `{cigar}`"
                )))
            }
        }
        length = 0;
        has_length = false;
    }
    if has_length {
        return Err(OrgraftError::InvalidArgument(format!(
            "trailing CIGAR length in `{cigar}`"
        )));
    }
    Ok((exact, mismatch))
}

fn parse_paf(line: &str) -> Result<PafHit, OrgraftError> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 12 {
        return Err(OrgraftError::InvalidArgument(format!(
            "malformed PAF line: {line}"
        )));
    }
    let matches = parse_usize_value(fields[9], "PAF matches")?;
    let block_length = parse_usize_value(fields[10], "PAF block length")?;
    let (substitution_matches, substitution_mismatches) = fields
        .iter()
        .find_map(|field| field.strip_prefix("cg:Z:"))
        .map(parse_eqx_counts)
        .transpose()?
        .unwrap_or_else(|| (matches, block_length.saturating_sub(matches)));
    let exact_flanks = fields
        .iter()
        .find_map(|field| field.strip_prefix("cg:Z:"))
        .and_then(|cigar| parse_cigar(cigar).ok())
        .is_some_and(|ops| {
            ops.first().is_some_and(|op| op.code == b'=' && op.count >= UNIQUE_PROJECTION_FLANK_BP)
                && ops.last().is_some_and(|op| op.code == b'=' && op.count >= UNIQUE_PROJECTION_FLANK_BP)
        });
    Ok(PafHit {
        query_name: fields[0].to_string(),
        query_length: parse_usize_value(fields[1], "PAF query length")?,
        query_start: parse_usize_value(fields[2], "PAF query start")?,
        query_end: parse_usize_value(fields[3], "PAF query end")?,
        strand: parse_strand(fields[4])?,
        target_start: parse_usize_value(fields[7], "PAF target start")?,
        target_end: parse_usize_value(fields[8], "PAF target end")?,
        matches,
        block_length,
        mapq: parse_usize_value(fields[11], "PAF mapq")?,
        substitution_matches,
        substitution_mismatches,
        primary: fields.iter().any(|field| *field == "tp:A:P"),
        exact_flanks,
    })
}

fn sort_hits(hits: &mut [PafHit]) {
    hits.sort_by(|a, b| {
        b.query_end
            .saturating_sub(b.query_start)
            .cmp(&a.query_end.saturating_sub(a.query_start))
            .then_with(|| b.matches.cmp(&a.matches))
            .then_with(|| {
                b.identity()
                    .partial_cmp(&a.identity())
                    .unwrap_or(Ordering::Equal)
            })
            .then_with(|| b.primary.cmp(&a.primary))
            .then_with(|| b.mapq.cmp(&a.mapq))
            .then_with(|| a.target_start.cmp(&b.target_start))
    });
}

fn dedupe_hits(mut hits: Vec<PafHit>, reference_len: usize) -> Vec<PafHit> {
    sort_hits(&mut hits);
    let mut out: Vec<PafHit> = Vec::new();
    for hit in hits {
        if out
            .iter()
            .any(|old| hits_share_circular_placement(old, &hit, reference_len))
        {
            continue;
        }
        out.push(hit);
    }
    out
}

fn hits_share_circular_placement(left: &PafHit, right: &PafHit, reference_len: usize) -> bool {
    if left.strand != right.strand {
        return false;
    }
    let query_overlap = left
        .query_end
        .min(right.query_end)
        .saturating_sub(left.query_start.max(right.query_start));
    let min_query_span = left
        .query_end
        .saturating_sub(left.query_start)
        .min(right.query_end.saturating_sub(right.query_start));
    if min_query_span == 0 || query_overlap.saturating_mul(100) < min_query_span.saturating_mul(90)
    {
        return false;
    }

    let left_folded = left.folded(reference_len);
    let right_folded = right.folded(reference_len);
    let target_overlap = folded_interval_overlap(&left_folded, &right_folded, reference_len);
    let min_target_span = folded_interval_len(&left_folded, reference_len)
        .min(folded_interval_len(&right_folded, reference_len));
    min_target_span > 0 && target_overlap.saturating_mul(100) >= min_target_span.saturating_mul(90)
}

fn folded_interval_overlap(
    left: &FoldedInterval,
    right: &FoldedInterval,
    reference_len: usize,
) -> usize {
    folded_interval_parts(left, reference_len)
        .into_iter()
        .flat_map(|left_part| {
            folded_interval_parts(right, reference_len)
                .into_iter()
                .map(move |right_part| {
                    interval_overlap(left_part.0, left_part.1, right_part.0, right_part.1)
                })
        })
        .sum()
}

fn folded_interval_len(interval: &FoldedInterval, reference_len: usize) -> usize {
    folded_interval_parts(interval, reference_len)
        .into_iter()
        .map(|(start, end)| end.saturating_sub(start) + 1)
        .sum()
}

fn folded_interval_parts(interval: &FoldedInterval, _reference_len: usize) -> Vec<(usize, usize)> {
    // Placement-overlap scores compare covered reference positions, not how many
    // times a long hit visits them. Export uses interval.parts without this union.
    merge_intervals(&interval.parts)
}

pub(crate) fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>, OrgraftError> {
    let text = fs::read_to_string(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
    })?;
    let mut records = Vec::new();
    let mut header: Option<String> = None;
    let mut sequence = String::new();
    for line in text.lines() {
        if let Some(next_header) = line.strip_prefix('>') {
            if let Some(old_header) = header.replace(next_header.trim().to_string()) {
                records.push(FastaRecord {
                    id: fasta_id(&old_header),
                    header: old_header,
                    sequence: sequence.clone(),
                });
                sequence.clear();
            }
        } else if !line.trim().is_empty() {
            sequence.push_str(&line.trim().to_ascii_uppercase());
        }
    }
    if let Some(old_header) = header {
        records.push(FastaRecord {
            id: fasta_id(&old_header),
            header: old_header,
            sequence,
        });
    }
    if records.is_empty() {
        return Err(OrgraftError::InvalidArgument(format!(
            "{} contains no FASTA records",
            path.display()
        )));
    }
    Ok(records)
}

fn write_single_fasta(path: &Path, header: &str, sequence: &str) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    write_fasta_record(&mut out, header, sequence)
}

fn write_node_fasta(path: &Path, gfa: &Gfa) -> Result<(), OrgraftError> {
    write_node_fasta_from_order(path, gfa)
}

fn write_node_fasta_from_order(path: &Path, gfa: &Gfa) -> Result<(), OrgraftError> {
    let mut out = File::create(path)?;
    for name in &gfa.order {
        let segment = &gfa.segments[name];
        write_fasta_record(&mut out, &segment.name, &segment.sequence)?;
    }
    Ok(())
}

pub(crate) fn write_fasta_record<W: Write>(
    out: &mut W,
    header: &str,
    sequence: &str,
) -> Result<(), OrgraftError> {
    writeln!(out, ">{header}")?;
    for chunk in sequence.as_bytes().chunks(80) {
        writeln!(out, "{}", String::from_utf8_lossy(chunk))?;
    }
    Ok(())
}

fn extract_circular(sequence: &str, start: usize, end: usize) -> String {
    let bytes = sequence.as_bytes();
    let len = bytes.len();
    let mut out = Vec::with_capacity(end.saturating_sub(start));
    for pos in start..end {
        out.push(bytes[pos % len]);
    }
    String::from_utf8_lossy(&out).to_string()
}

fn normalize_repeat_constrained_sequence(
    candidate: &str,
    reference: &str,
    boundary_offsets: &[usize],
) -> Result<(String, char, usize, usize), OrgraftError> {
    if candidate.len() != reference.len() {
        return Err(OrgraftError::InvalidArgument(format!(
            "repeat-constrained candidate length {} differs from polished input length {}; validate-only handoff requires equal lengths",
            candidate.len(),
            reference.len()
        )));
    }
    if candidate.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "cannot normalize an empty repeat-constrained sequence".to_string(),
        ));
    }

    let reverse_candidate = reverse_complement(candidate);
    let mut plus_offsets = boundary_offsets
        .iter()
        .map(|offset| offset % candidate.len())
        .collect::<Vec<_>>();
    plus_offsets.push(0);
    plus_offsets.sort_unstable();
    plus_offsets.dedup();
    let mut minus_offsets = plus_offsets
        .iter()
        .map(|offset| (candidate.len() - offset) % candidate.len())
        .collect::<Vec<_>>();
    minus_offsets.sort_unstable();
    minus_offsets.dedup();

    let mut best: Option<(usize, char, usize, String)> = None;
    for (orientation, oriented_candidate) in [('+', candidate), ('-', reverse_candidate.as_str())] {
        let offsets = if orientation == '+' {
            &plus_offsets
        } else {
            &minus_offsets
        };
        for rotation in offsets {
            let normalized = rotate_sequence_left(oriented_candidate, *rotation);
            let exact_matches = normalized
                .bytes()
                .zip(reference.bytes())
                .filter(|(left, right)| left == right)
                .count();
            let candidate_key = (exact_matches, orientation, *rotation, normalized);
            let replace = best.as_ref().is_none_or(|current| {
                candidate_key.0 > current.0
                    || (candidate_key.0 == current.0 && candidate_key.1 == '+' && current.1 == '-')
                    || (candidate_key.0 == current.0
                        && candidate_key.1 == current.1
                        && candidate_key.2 < current.2)
            });
            if replace {
                best = Some(candidate_key);
            }
        }
    }

    best.map(|(matches, orientation, rotation, normalized)| {
        (normalized, orientation, rotation, matches)
    })
    .ok_or_else(|| {
        OrgraftError::InvalidArgument(
            "repeat-constrained candidate has no graph-node boundary available for circular normalization"
                .to_string(),
        )
    })
}

fn rotate_sequence_left(sequence: &str, offset: usize) -> String {
    let offset = offset % sequence.len();
    format!("{}{}", &sequence[offset..], &sequence[..offset])
}

fn reverse_complement(sequence: &str) -> String {
    sequence
        .bytes()
        .rev()
        .map(|base| match base {
            b'A' | b'a' => 'T',
            b'C' | b'c' => 'G',
            b'G' | b'g' => 'C',
            b'T' | b't' => 'A',
            b'N' | b'n' => 'N',
            other => other as char,
        })
        .collect()
}

fn complement_base(base: char) -> char {
    match base {
        'A' | 'a' => 'T',
        'C' | 'c' => 'G',
        'G' | 'g' => 'C',
        'T' | 't' => 'A',
        'N' | 'n' => 'N',
        other => other,
    }
}

fn flip_orient(orient: char) -> char {
    if orient == '+' {
        '-'
    } else {
        '+'
    }
}

fn parse_orient(value: &str, path: &Path, line: usize) -> Result<char, OrgraftError> {
    match value {
        "+" => Ok('+'),
        "-" => Ok('-'),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "{}:{} invalid orientation `{value}`",
            path.display(),
            line
        ))),
    }
}

fn parse_strand(value: &str) -> Result<char, OrgraftError> {
    match value {
        "+" => Ok('+'),
        "-" => Ok('-'),
        _ => Err(OrgraftError::InvalidArgument(format!(
            "invalid strand `{value}`"
        ))),
    }
}

pub(crate) fn required_value<'a>(
    args: &'a [String],
    index: &mut usize,
    option: &str,
) -> Result<&'a str, OrgraftError> {
    *index += 1;
    args.get(*index)
        .map(String::as_str)
        .ok_or_else(|| OrgraftError::InvalidArgument(format!("{option} requires a value")))
}

pub(crate) fn parse_usize(value: &str, option: &str) -> Result<usize, OrgraftError> {
    value.parse::<usize>().map_err(|_| {
        OrgraftError::InvalidArgument(format!("{option} must be an integer, got `{value}`"))
    })
}

fn parse_subgraph_id(value: &str) -> Result<String, OrgraftError> {
    if value.is_empty() {
        return Err(OrgraftError::InvalidArgument(
            "--subgraph must not be empty".to_string(),
        ));
    }
    if !value
        .chars()
        .all(|ch| ch.is_ascii_alphanumeric() || ch == '_' || ch == '-')
    {
        return Err(OrgraftError::InvalidArgument(format!(
            "--subgraph may only contain ASCII letters, digits, `_`, and `-`, got `{value}`"
        )));
    }
    Ok(value.to_string())
}

fn parse_usize_value(value: &str, label: &str) -> Result<usize, OrgraftError> {
    value.parse::<usize>().map_err(|_| {
        OrgraftError::InvalidArgument(format!("{label} must be an integer, got `{value}`"))
    })
}

pub(crate) fn canonicalize_existing(path: &Path, label: &str) -> Result<PathBuf, OrgraftError> {
    fs::canonicalize(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {label} {}: {error}", path.display()))
    })
}

fn remove_path_if_exists(path: &Path) -> Result<(), OrgraftError> {
    match fs::symlink_metadata(path) {
        Ok(metadata) if metadata.is_dir() => fs::remove_dir_all(path).map_err(OrgraftError::Io),
        Ok(_) => fs::remove_file(path).map_err(OrgraftError::Io),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(OrgraftError::Io(error)),
    }
}

fn read_soft_paths(path: &Path) -> Result<HashMap<String, PathBuf>, OrgraftError> {
    if !path.exists() {
        return Ok(HashMap::new());
    }
    let text = fs::read_to_string(path).map_err(|error| {
        OrgraftError::InvalidArgument(format!("cannot read {}: {error}", path.display()))
    })?;
    let mut map = HashMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 2 {
            map.insert(fields[0].to_string(), PathBuf::from(fields[1]));
        }
    }
    Ok(map)
}

fn discover_validate_data_dir(polished_fasta: &Path) -> Option<PathBuf> {
    let file_dir = polished_fasta.parent()?;
    let round_dir = match file_dir.file_name().and_then(|name| name.to_str()) {
        Some("01.inputs") | Some("02.polish") => file_dir.parent()?,
        _ => return None,
    };
    let candidate = round_dir.join("03.validate/01.data");
    if candidate.is_dir() {
        Some(candidate)
    } else {
        None
    }
}

fn fasta_id(header: &str) -> String {
    header
        .split_whitespace()
        .next()
        .unwrap_or(header)
        .to_string()
}

fn replace_tag(tags: &[String], new_tag: String) -> Vec<String> {
    let tag_name = new_tag.split(':').next().unwrap_or("");
    let prefix = format!("{tag_name}:");
    let mut out: Vec<String> = tags
        .iter()
        .filter(|tag| !tag.starts_with(&prefix))
        .cloned()
        .collect();
    out.push(new_tag);
    out
}

fn replace_tags(tags: &[String], new_tags: &[String]) -> Vec<String> {
    let mut out = tags.to_vec();
    for tag in new_tags {
        out = replace_tag(&out, tag.clone());
    }
    out
}

fn tag_value(tags: &[String], tag_name: &str) -> Option<String> {
    let prefix = format!("{tag_name}:");
    tags.iter().find_map(|tag| {
        if !tag.starts_with(&prefix) {
            return None;
        }
        tag.splitn(3, ':').nth(2).map(str::to_string)
    })
}

fn numeric_tag(tags: &[String], tag_name: &str) -> Option<f64> {
    tag_value(tags, tag_name).and_then(|value| value.parse::<f64>().ok())
}

fn link_support(link: &Link) -> f64 {
    for tag_name in ["RC", "SK", "PA"] {
        if let Some(value) = numeric_tag(&link.tags, tag_name) {
            return value;
        }
    }
    0.0
}

// pub(crate): physical-end side of a node on a link (refine uses the
// identical rule to recognise an unresolved 2-left x 2-right graph repeat.
pub(crate) fn link_endpoint_side(link: &Link, node: &str) -> Option<char> {
    if link.from_name == node {
        Some(if link.from_orient == '-' { '-' } else { '+' })
    } else if link.to_name == node {
        Some(if link.to_orient == '-' { '+' } else { '-' })
    } else {
        None
    }
}

fn natural_cmp(left: &str, right: &str) -> Ordering {
    natural_key(left).cmp(&natural_key(right))
}

fn natural_key(value: &str) -> Vec<NaturalPart> {
    let mut parts = Vec::new();
    let mut current = String::new();
    let mut in_digits = None;
    for ch in value.chars() {
        let digit = ch.is_ascii_digit();
        match in_digits {
            Some(state) if state == digit => current.push(ch),
            Some(state) => {
                parts.push(if state {
                    NaturalPart::Number(current.parse::<u64>().unwrap_or(0))
                } else {
                    NaturalPart::Text(current.clone())
                });
                current.clear();
                current.push(ch);
                in_digits = Some(digit);
            }
            None => {
                current.push(ch);
                in_digits = Some(digit);
            }
        }
    }
    if let Some(state) = in_digits {
        parts.push(if state {
            NaturalPart::Number(current.parse::<u64>().unwrap_or(0))
        } else {
            NaturalPart::Text(current)
        });
    }
    parts
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum NaturalPart {
    Text(String),
    Number(u64),
}

fn stable_hash64(sequence: &str) -> String {
    let mut hash: u64 = 0xcbf29ce484222325;
    for byte in sequence.as_bytes() {
        hash ^= *byte as u64;
        hash = hash.wrapping_mul(0x100000001b3);
    }
    format!("{hash:016x}")
}

fn temp_work_dir(label: &str) -> PathBuf {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|value| value.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!("{label}-{}-{stamp}", process::id()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recorded_path_options_require_paired_evidence_and_fixed_topology() {
        let base = [
            "--edited-gfa",
            "component.gfa",
            "--polished-fasta",
            "final.fasta",
        ];
        let parse = |extra: &[&str]| {
            RebuildOptions::from_args(
                &base
                    .iter()
                    .chain(extra)
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>(),
            )
        };
        assert!(parse(&[]).is_ok());
        assert!(parse(&["--checked-gfa", "checked.gfa"]).is_err());
        assert!(parse(&["--resolve-details", "details.tsv"]).is_err());
        let paired = [
            "--checked-gfa",
            "checked.gfa",
            "--resolve-details",
            "details.tsv",
        ];
        assert!(parse(&paired).is_ok());
        assert!(parse(&["--ordinary-source-dir", "round_3"]).is_err());
        let mut ordinary = paired.to_vec();
        ordinary.extend(["--ordinary-source-dir", "round_3"]);
        assert!(parse(&ordinary).is_err());
        ordinary.extend(["--reads", "reads.fastq.gz"]);
        assert!(parse(&ordinary).is_ok());
        for flag in ["--force", "--overlap-only"] {
            let mut conflict = ordinary.clone(); conflict.push(flag);
            assert!(parse(&conflict).is_err());
        }
        for conflict in ["--repeat-resolution-node", "--merged-gfa-template"] {
            let mut args = paired.to_vec();
            args.extend([conflict, "value"]);
            assert!(parse(&args).is_err());
        }
    }

    #[test]
    fn parses_options() {
        let args = vec![
            "--edited-gfa".to_string(),
            "graph.gfa".to_string(),
            "--polished-fasta".to_string(),
            "final.fa".to_string(),
            "--reads".to_string(),
            "reads.fastq.gz".to_string(),
            "--image-reference-fasta".to_string(),
            "bait.fa".to_string(),
            "--out-dir".to_string(),
            "out".to_string(),
        ];
        let options = RebuildOptions::from_args(&args).unwrap();
        assert_eq!(options.edited_gfa, PathBuf::from("graph.gfa"));
        assert_eq!(options.polished_fasta, PathBuf::from("final.fa"));
        assert_eq!(options.reads, Some(PathBuf::from("reads.fastq.gz")));
        assert_eq!(
            options.image_reference_fasta,
            Some(PathBuf::from("bait.fa"))
        );
        assert_eq!(options.out_dir, PathBuf::from("out"));
        assert_eq!(options.subgraph, "subgraph_001");
    }

    #[test]
    fn rebuild_out_dir_defaults_to_results_rebuild() {
        let args = vec![
            "--edited-gfa".to_string(),
            "graph.gfa".to_string(),
            "--polished-fasta".to_string(),
            "final.fa".to_string(),
        ];
        let options = RebuildOptions::from_args(&args).unwrap();
        assert_eq!(options.out_dir, PathBuf::from("results/rebuild"));
        assert_eq!(options.image_reference_fasta, None);
    }

    #[test]
    fn image_reference_fasta_is_optional() {
        let args = vec![
            "--edited-gfa".to_string(),
            "graph.gfa".to_string(),
            "--polished-fasta".to_string(),
            "final.fa".to_string(),
        ];
        let options = RebuildOptions::from_args(&args).unwrap();
        assert_eq!(options.image_reference_fasta, None);
    }

    #[test]
    fn discover_validate_data_dir_uses_workflow_round_directory() {
        let base = std::env::temp_dir().join(format!(
            "orgraft-rebuild-round-suffix-{}",
            std::process::id()
        ));
        let round_dir = base.join("04.polish/mito/subgraph_001/round_2");
        let data_dir = round_dir.join("03.validate/01.data");
        fs::create_dir_all(&data_dir).unwrap();
        let polished = round_dir.join("01.inputs/linear_subgraph.round_2.fasta");
        fs::create_dir_all(polished.parent().unwrap()).unwrap();
        fs::write(&polished, ">x\nACGT\n").unwrap();

        assert_eq!(discover_validate_data_dir(&polished), Some(data_dir));
        let _ = fs::remove_dir_all(base);
    }

    #[test]
    fn gfa_image_export_failure_is_recorded_without_error() {
        let row = run_gfa_editor_image(
            Path::new("/definitely/missing/gfa_editor_cli"),
            &HashMap::new(),
            Path::new("input.gfa"),
            Path::new("output.pdf"),
            Path::new("reference.fa"),
            "pdf",
        );

        assert_eq!(row.format, "pdf");
        assert_eq!(row.status, "failed_to_run");
        assert!(row.stderr.contains("failed to run GFA_Editor image export"));
    }

    #[test]
    fn skipped_gfa_reference_images_records_both_formats() {
        let paths = OutputPaths::new(Path::new("out"), "subgraph_001");
        let rows = skipped_gfa_reference_images(&paths, "missing cli");

        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].format, "pdf");
        assert_eq!(rows[0].status, "skipped_missing_gfa_editor_cli");
        assert_eq!(rows[0].stderr, "missing cli");
        assert_eq!(rows[1].format, "svg");
        assert_eq!(rows[1].status, "skipped_missing_gfa_editor_cli");
    }

    #[test]
    fn help_places_subgraph_in_inputs_and_uses_layout_summary() {
        assert!(HELP.contains(
            "Usage:\n  orgraft rebuild --organelle NAME --subgraph ID --edited-gfa FILE --polished-fasta FILE [options]"
        ));
        assert!(!HELP.contains(
            "Usage:\n  orgraft rebuild --organelle NAME --subgraph ID --edited-gfa FILE --polished-fasta FILE --out-dir DIR"
        ));
        assert!(HELP.contains(
            "  --out-dir DIR                 rebuild root output directory [results/rebuild]\n  --force                       replace only this component's outputs"
        ));
        assert!(HELP.contains(
            "Inputs:\n  --organelle NAME              organelle name for this rebuild run [mito]\n  --subgraph ID                 subgraph/ring id [subgraph_001]"
        ));
        assert!(HELP.contains(
            "Additional Parameters:\n  --threads N                   minimap2 threads for node-to-polished projection [4]\n  --image-reference-fasta FILE  reference FASTA for graph colouring; enables PDF/SVG export"
        ));
        assert!(HELP.contains(
            "Layout: OUT/SUBGRAPH/rebuild_SUBGRAPH* plus OUT/logs/*.tsv (pdf/svg need reference)"
        ));
        assert!(!HELP.contains("DIR/logs/manifest.tsv"));
        assert!(!HELP.contains("--genome"));
    }

    #[test]
    fn extracts_circular_sequence() {
        assert_eq!(extract_circular("ABCDEF", 4, 8), "EFAB");
    }

    #[test]
    fn repeat_path_endpoints_are_oriented_for_repeat_plus_traversal() {
        let incoming_reverse = Link {
            from_name: "repeat".to_string(),
            from_orient: '-',
            to_name: "left".to_string(),
            to_orient: '+',
            overlap: "0M".to_string(),
            tags: Vec::new(),
        };
        let outgoing_reverse = Link {
            from_name: "right".to_string(),
            from_orient: '+',
            to_name: "repeat".to_string(),
            to_orient: '-',
            overlap: "0M".to_string(),
            tags: Vec::new(),
        };

        assert_eq!(link_endpoint_side(&incoming_reverse, "repeat"), Some('-'));
        assert_eq!(
            repeat_path_edge(&incoming_reverse, "repeat", 1)
                .unwrap()
                .label(),
            "left-"
        );
        assert_eq!(link_endpoint_side(&outgoing_reverse, "repeat"), Some('+'));
        assert_eq!(
            repeat_path_edge(&outgoing_reverse, "repeat", 2)
                .unwrap()
                .label(),
            "right-"
        );
    }

    #[test]
    fn repeat_path_records_use_integer_spanning_read_support() {
        let mut gfa = test_gfa(&[("left", "AC"), ("repeat", "GT"), ("right", "TA")]);
        gfa.segments.get_mut("repeat").unwrap().tags =
            vec!["PM:Z:repeat_path_support_pending".to_string()];
        let rows = vec![RepeatPathSupportRow {
            repeat_node: "repeat".to_string(),
            repeat_status: "unresolved".to_string(),
            support_method: "spanning_reads".to_string(),
            left_edge_index: Some(1),
            right_edge_index: Some(2),
            left_endpoint: "left+".to_string(),
            right_endpoint: "right+".to_string(),
            left_physical_endpoint: "left+".to_string(),
            right_physical_endpoint: "right+".to_string(),
            left_support: Some(100.0),
            right_support: Some(90.0),
            path_support: Some(17),
            path_ratio: Some(1.0),
            supporting_read_ids: vec!["read_1".to_string()],
        }];

        attach_repeat_path_support_paths(&mut gfa, &rows);

        assert!(gfa.segments["repeat"]
            .tags
            .contains(&"PM:Z:spanning_reads".to_string()));
        assert_eq!(gfa.other_lines.len(), 1);
        assert!(gfa.other_lines[0].contains(&"RC:i:17".to_string()));
        assert!(gfa.other_lines[0].contains(&"PM:Z:spanning_reads".to_string()));
    }

    #[test]
    fn dominant_repeat_pairing_uses_complete_spanning_read_counts() {
        let make_row = |left_edge_index,
                        right_edge_index,
                        left_endpoint: &str,
                        right_endpoint: &str,
                        path_support| RepeatPathSupportRow {
            repeat_node: "utg5".to_string(),
            repeat_status: "unresolved".to_string(),
            support_method: "spanning_reads".to_string(),
            left_edge_index: Some(left_edge_index),
            right_edge_index: Some(right_edge_index),
            left_endpoint: left_endpoint.to_string(),
            right_endpoint: right_endpoint.to_string(),
            left_physical_endpoint: left_endpoint.to_string(),
            right_physical_endpoint: right_endpoint.to_string(),
            left_support: Some(path_support as f64),
            right_support: Some(path_support as f64),
            path_support: Some(path_support),
            path_ratio: None,
            supporting_read_ids: Vec::new(),
        };
        let rows = vec![
            make_row(1, 3, "long+", "utg1+", 242),
            make_row(1, 4, "long+", "utg7+", 4),
            make_row(2, 3, "utg2+", "utg1+", 6),
            make_row(2, 4, "utg2+", "utg7+", 223),
        ];

        let decision = select_dominant_repeat_pairing(&rows, "utg5").unwrap();

        assert_eq!(decision.selected_paths, vec!["p1", "p4"]);
        assert_eq!(decision.selected_support, 465);
        assert_eq!(decision.alternative_support, 10);
        assert_eq!(decision.pairing_signature, "long+--utg1+|utg2+--utg7+");
    }

    #[test]
    fn final_rebuild_preserves_the_initial_candidate_id_as_provenance() {
        let base = temp_work_dir("orgraft-repeat-candidate-provenance");
        let paths = OutputPaths::new(&base, "subgraph_001");
        paths.create(false).unwrap();
        fs::write(
            &paths.result_stats,
            "section\tsubgraph\titem\tmetric\tvalue\textra\nsummary\tsubgraph_001\t.\trepeat_resolution_node\tutg5\t.\nsummary\tsubgraph_001\t.\trepeat_resolution_candidate\tauto_repeat_003\t.\n",
        )
        .unwrap();
        fs::write(
            &paths.repeat_resolution,
            "section\titem\tpairing_signature\tselected\tscore\tcontinuous_fraction\tdetails\nsummary\tutg5\tlong+--utg1+|utg2+--utg7+\tYes\t.\t.\tmethod=unique_after_pairing_filter\n",
        )
        .unwrap();

        let prior = read_prior_repeat_resolution(&paths, "utg5").unwrap();
        let decision = RepeatResolutionDecision {
            repeat_node: "utg5".to_string(),
            pairing_signature: "long+--utg1+|utg2+--utg7+".to_string(),
            selected_paths: vec!["p1".to_string(), "p4".to_string()],
            selected_support: 482,
            alternative_support: 12,
        };

        assert_eq!(
            source_repeat_candidate(Some(&prior), &decision, "auto_repeat_002"),
            "auto_repeat_003"
        );
        let mut changed = decision.clone();
        changed.pairing_signature = "long+--utg7+|utg2+--utg1+".to_string();
        assert_eq!(
            source_repeat_candidate(Some(&prior), &changed, "auto_repeat_002"),
            "auto_repeat_002"
        );

        let _ = fs::remove_dir_all(base);
    }

    #[test]
    fn repeat_path_identity_does_not_penalize_template_only_gap_length() {
        let alignment = RepeatPathAlignment {
            read_id: "read_1".to_string(),
            target_id: "repeat_path_0001".to_string(),
            query_start: 3_436,
            query_end: 7_876,
            target_start: 0,
            target_end: 4_953,
            matches: 4_436,
            mapq: 60,
            primary: true,
        };

        assert!(alignment.identity() > 0.99);
        assert!(alignment.identity() >= REPEAT_PATH_MIN_IDENTITY);
    }

    #[test]
    fn repeat_candidate_normalization_keeps_a_graph_node_boundary() {
        let candidate = "AAAACCGT";
        let reference = "AACCGTAA";

        let (normalized, _, rotation, exact_matches) =
            normalize_repeat_constrained_sequence(candidate, reference, &[0, 4]).unwrap();

        assert!(matches!(rotation, 0 | 4));
        assert_ne!(normalized, reference);
        assert!(exact_matches < candidate.len());
    }

    #[test]
    fn repeat_core_alignment_threshold_adapts_below_one_kilobase() {
        assert_eq!(repeat_core_alignment_threshold(873), 655);
        assert_eq!(repeat_core_alignment_threshold(2_000), 1_000);
        assert_eq!(repeat_core_alignment_threshold(400), 500);
    }

    #[test]
    fn repeat_core_assignment_requires_two_mapping_anchored_copies() {
        let pair = RepeatCorePair {
            length: 878,
            source_blast_length: 878,
            source_pident: 99.886,
            source_gapopen: 1,
            source_mismatch: 0,
            boundary_extension_left: 0,
            boundary_extension_right: 0,
            copies: vec![
                RepeatCopy {
                    start: 1_000,
                    end: 1_877,
                    strand: '+',
                },
                RepeatCopy {
                    start: 5_000,
                    end: 5_877,
                    strand: '-',
                },
            ],
        };
        let both_copies = test_mapping(vec![
            make_synthetic_hit("repeat", 878, 1, 878, 1_000, 1_877, '+'),
            make_synthetic_hit("repeat", 878, 1, 878, 5_000, 5_877, '-'),
        ]);
        let doubled_reference_copies = test_mapping(vec![
            make_synthetic_hit("repeat", 878, 1, 878, 21_000, 21_877, '+'),
            make_synthetic_hit("repeat", 878, 1, 878, 25_000, 25_877, '-'),
        ]);
        let one_copy = test_mapping(vec![make_synthetic_hit(
            "repeat", 878, 1, 878, 1_000, 1_877, '+',
        )]);
        let decoy = RepeatCorePair {
            copies: vec![
                RepeatCopy {
                    start: 10_000,
                    end: 10_877,
                    strand: '+',
                },
                RepeatCopy {
                    start: 15_000,
                    end: 15_877,
                    strand: '-',
                },
            ],
            ..pair.clone()
        };

        assert_eq!(
            repeat_pair_mapping_score(&both_copies, &pair, 20_000),
            Some(1_756)
        );
        assert_eq!(
            repeat_pair_mapping_score(&doubled_reference_copies, &pair, 20_000),
            Some(1_756)
        );
        assert_eq!(repeat_pair_mapping_score(&one_copy, &pair, 20_000), None);
        assert_eq!(
            repeat_pair_mapping_score(&both_copies, &decoy, 20_000),
            None
        );
    }

    #[test]
    fn repeat_core_assignment_rescues_indel_heavy_two_copy_anchors() {
        let pair = RepeatCorePair {
            length: 1_000,
            source_blast_length: 1_000,
            source_pident: 100.0,
            source_gapopen: 0,
            source_mismatch: 0,
            boundary_extension_left: 0,
            boundary_extension_right: 0,
            copies: vec![
                RepeatCopy {
                    start: 1_000,
                    end: 1_999,
                    strand: '+',
                },
                RepeatCopy {
                    start: 5_000,
                    end: 5_999,
                    strand: '-',
                },
            ],
        };
        let make_indel_heavy_hit = |target_start, target_end, strand| PafHit {
            query_name: "repeat".to_string(),
            query_length: 1_000,
            query_start: 0,
            query_end: 1_000,
            strand,
            target_start,
            target_end,
            matches: 979,
            block_length: 1_000,
            mapq: 0,
            substitution_matches: 979,
            substitution_mismatches: 1,
            primary: true,
            exact_flanks: false,
        };
        let hits = vec![
            make_indel_heavy_hit(999, 1_979, '+'),
            make_indel_heavy_hit(4_999, 5_979, '-'),
        ];
        let mapping = NodeMapping {
            selected: hits.first().cloned(),
            accepted_hits: Vec::new(),
            all_hits: hits,
        };

        assert_eq!(
            repeat_pair_mapping_score(&mapping, &pair, 20_000),
            Some(1_960)
        );
    }

    #[test]
    fn repeat_core_assignment_rejects_substitution_heavy_relaxed_hits() {
        let pair = RepeatCorePair {
            length: 1_000,
            source_blast_length: 1_000,
            source_pident: 100.0,
            source_gapopen: 0,
            source_mismatch: 0,
            boundary_extension_left: 0,
            boundary_extension_right: 0,
            copies: vec![
                RepeatCopy {
                    start: 1_000,
                    end: 1_999,
                    strand: '+',
                },
                RepeatCopy {
                    start: 5_000,
                    end: 5_999,
                    strand: '-',
                },
            ],
        };
        let make_substitution_heavy_hit = |target_start, target_end, strand| PafHit {
            query_name: "repeat".to_string(),
            query_length: 1_000,
            query_start: 0,
            query_end: 1_000,
            strand,
            target_start,
            target_end,
            matches: 979,
            block_length: 1_000,
            mapq: 0,
            substitution_matches: 979,
            substitution_mismatches: 21,
            primary: true,
            exact_flanks: false,
        };
        let hits = vec![
            make_substitution_heavy_hit(999, 1_979, '+'),
            make_substitution_heavy_hit(4_999, 5_979, '-'),
        ];
        let mapping = NodeMapping {
            selected: hits.first().cloned(),
            accepted_hits: Vec::new(),
            all_hits: hits,
        };

        assert_eq!(repeat_pair_mapping_score(&mapping, &pair, 20_000), None);
    }

    #[test]
    fn circular_projection_dedupes_wrap_and_truncated_same_placement() {
        let make_hit = |query_start, query_end, target_start, target_end| PafHit {
            query_name: "long_node".to_string(),
            query_length: 243_451,
            query_start,
            query_end,
            strand: '+',
            target_start,
            target_end,
            matches: query_end - query_start,
            block_length: query_end - query_start,
            mapq: 60,
            substitution_matches: query_end - query_start,
            substitution_mismatches: 0,
            primary: query_start == 0,
            exact_flanks: false,
        };
        let reference_len = 339_849;
        let wrap = make_hit(0, 243_451, 339_554, 580_840);
        let truncated = make_hit(295, 243_451, 0, 240_991);

        assert!(hits_share_circular_placement(
            &wrap,
            &truncated,
            reference_len
        ));
        assert_eq!(dedupe_hits(vec![truncated, wrap], reference_len).len(), 1);

        let first_copy = make_hit(0, 1_000, 10_000, 11_000);
        let second_copy = make_hit(0, 1_000, 40_000, 41_000);
        assert!(!hits_share_circular_placement(
            &first_copy,
            &second_copy,
            reference_len
        ));
    }

    #[test]
    fn folded_intervals_keep_short_spans_complete_turns_and_extra_turns() {
        for (start, end, expected) in [
            (2, 6, "3-6"),
            (8, 12, "9-10,1-2"),
            (0, 10, "1-10"),
            (3, 13, "4-10,1-3"),
            (8, 22, "9-10,1-10,1-2"),
            (23, 27, "4-7"),
        ] {
            let interval = FoldedInterval::new(start, end, 10);
            assert_eq!(interval.parts_text(10), expected);
            assert_eq!(
                interval.parts.iter().map(|(a, b)| b - a + 1).sum::<usize>(),
                end - start
            );
        }
        let long = FoldedInterval::new(8, 22, 10);
        let full = FoldedInterval::new(0, 10, 10);
        assert_eq!(folded_interval_len(&long, 10), 10);
        assert_eq!(folded_interval_overlap(&long, &full, 10), 10);
    }

    #[test]
    fn folded_remapped_depth_keeps_crossing_positions_without_double_counting_turns() {
        let base = temp_work_dir("orgraft-circular-depth-test");
        fs::create_dir_all(&base).unwrap();
        let mut depth = "position\tfl_depth\n".to_string();
        for position in 1..=10 {
            writeln!(depth, "{position}\t{position}").unwrap();
        }
        fs::write(base.join("sv_coverage.tsv"), depth).unwrap();
        let mapping = Mapping {
            by_node: BTreeMap::from([
                (
                    "crossing".to_string(),
                    test_mapping(vec![make_synthetic_hit("crossing", 4, 1, 4, 9, 12, '+')]),
                ),
                (
                    "extra_turn".to_string(),
                    test_mapping(vec![make_synthetic_hit(
                        "extra_turn",
                        14,
                        1,
                        14,
                        9,
                        22,
                        '+',
                    )]),
                ),
            ]),
            reference_len: 10,
            command: String::new(),
            stderr: String::new(),
        };
        let classes = BTreeMap::from([
            ("crossing".to_string(), "single_copy_node".to_string()),
            ("extra_turn".to_string(), "single_copy_node".to_string()),
        ]);
        let rows = compute_node_remapped_coverage(&mapping, &classes, Some(&base)).unwrap();
        let crossing = rows.iter().find(|row| row.node == "crossing").unwrap();
        assert_eq!(crossing.intervals, "1-2,9-10");
        assert_eq!(crossing.covered_bases, Some(4));
        assert_eq!(crossing.depth_bases, Some(22.0));
        assert_eq!(crossing.depth_mean, Some(5.5));
        let extra = rows.iter().find(|row| row.node == "extra_turn").unwrap();
        assert_eq!(extra.intervals, "1-10");
        assert_eq!(extra.covered_bases, Some(10));
        assert_eq!(extra.depth_bases, Some(55.0));
        assert_eq!(
            mapping.by_node["extra_turn"].accepted_hits[0].target_end,
            22
        );
        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn folded_long_node_exports_every_turn_in_gfa_and_extract() {
        let base = temp_work_dir("orgraft-circular-extract-test");
        fs::create_dir_all(&base).unwrap();
        let sequence = "TAACGTACGTAACG";
        let gfa = test_gfa(&[("node", sequence)]);
        let mapping = Mapping {
            by_node: BTreeMap::from([(
                "node".to_string(),
                test_mapping(vec![make_synthetic_hit("node", 14, 1, 14, 9, 22, '+')]),
            )]),
            reference_len: 10,
            command: String::new(),
            stderr: String::new(),
        };
        let classes = BTreeMap::from([("node".to_string(), "single_copy_node".to_string())]);
        let sequences = BTreeMap::from([("node".to_string(), sequence.to_string())]);
        let verified =
            make_verified_gfa_from_sequences(&gfa, &sequences, &mapping, &classes, "test");
        assert!(verified.segments["node"]
            .tags
            .contains(&"VS:Z:9-10,1-10,1-2:+".to_string()));
        let path = base.join("extract.tsv");
        write_rebuild_extract(
            &path,
            "subgraph_001",
            &verified,
            &mapping,
            &classes,
            &[],
            &[],
        )
        .unwrap();
        let text = fs::read_to_string(path).unwrap();
        let fields = text
            .lines()
            .skip(1)
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .collect::<Vec<_>>();
        assert_eq!(fields.len(), 3);
        assert_eq!(
            fields
                .iter()
                .map(|row| row[3].parse::<usize>().unwrap() - row[2].parse::<usize>().unwrap() + 1)
                .sum::<usize>(),
            14
        );
        assert!(fields
            .iter()
            .all(|row| row[10] == "1" && row.last().unwrap().contains("part=")));
        fs::remove_dir_all(base).unwrap();
    }

    #[test]
    fn consistency_uses_all_accepted_node_copies_for_linear_tiling() {
        let record = test_record("ACGTACGTAA");
        let gfa = test_gfa(&[("repeat", "AC"), ("single", "GTGTAA")]);
        let mut by_node = BTreeMap::new();
        let repeat_hits = vec![
            make_synthetic_hit("repeat", 2, 1, 2, 1, 2, '+'),
            make_synthetic_hit("repeat", 2, 1, 2, 5, 6, '+'),
        ];
        let single_hits = vec![
            make_synthetic_hit("single", 6, 1, 2, 3, 4, '+'),
            make_synthetic_hit("single", 6, 3, 6, 7, 10, '+'),
        ];
        by_node.insert("repeat".to_string(), test_mapping(repeat_hits));
        by_node.insert("single".to_string(), test_mapping(single_hits));
        let mapping = Mapping {
            by_node,
            reference_len: record.sequence.len(),
            command: String::new(),
            stderr: String::new(),
        };

        let consistency = CoordinateConsistency::new(&record, &gfa, &gfa, &gfa, &mapping);

        assert_eq!(consistency.covered_bases, 10);
        assert_eq!(consistency.gap_bases, 0);
        assert_eq!(consistency.multi_covered_bases, 0);
        assert_eq!(consistency.linear_tiling_status(), "PASS");
        assert_eq!(consistency.coverage_fraction, 1.0);
    }

    #[test]
    fn consistency_reports_gaps_and_overlaps() {
        let record = test_record("ACGTAC");
        let gfa = test_gfa(&[("left", "ACG"), ("right", "GT")]);
        let mut by_node = BTreeMap::new();
        by_node.insert(
            "left".to_string(),
            test_mapping(vec![make_synthetic_hit("left", 3, 1, 3, 1, 3, '+')]),
        );
        by_node.insert(
            "right".to_string(),
            test_mapping(vec![make_synthetic_hit("right", 2, 1, 2, 3, 4, '+')]),
        );
        let mapping = Mapping {
            by_node,
            reference_len: record.sequence.len(),
            command: String::new(),
            stderr: String::new(),
        };

        let consistency = CoordinateConsistency::new(&record, &gfa, &gfa, &gfa, &mapping);

        assert_eq!(consistency.covered_bases, 4);
        assert_eq!(consistency.gap_bases, 2);
        assert_eq!(consistency.multi_covered_bases, 1);
        assert_eq!(consistency.linear_tiling_status(), "WARN");
    }

    #[test]
    fn premerge_preserves_one_bidirected_self_link() {
        let mut raw = test_gfa(&[("u0", "ACGT")]);
        raw.links = vec![
            Link {
                from_name: "u0".to_string(),
                from_orient: '+',
                to_name: "u0".to_string(),
                to_orient: '+',
                overlap: "3M".to_string(),
                tags: vec!["EC:i:499".to_string()],
            },
            Link {
                from_name: "u0".to_string(),
                from_orient: '-',
                to_name: "u0".to_string(),
                to_orient: '-',
                overlap: "3M".to_string(),
                tags: vec!["EC:i:499".to_string()],
            },
        ];

        let (merged, mode) = merge_unambiguous_gfa(&raw);

        assert_eq!(mode, "input_already_merged_or_no_linear_compaction");
        assert_eq!(merged.order, vec!["u0"]);
        assert_eq!(merged.links.len(), 1);
        assert_eq!(merged.links[0].from_name, "u0");
        assert_eq!(merged.links[0].to_name, "u0");
        assert_eq!(merged.links[0].overlap, "3M");
        assert_eq!(merged.links[0].tags, vec!["EC:i:499"]);
        assert_eq!(merged.degrees()["u0"], 2);
        assert!(repeat_like_nodes(&merged).is_empty());
    }

    #[test]
    fn premerge_still_drops_links_absorbed_by_linear_compaction() {
        let mut raw = test_gfa(&[("a", "AC"), ("b", "GT")]);
        raw.links.push(Link {
            from_name: "a".to_string(),
            from_orient: '+',
            to_name: "b".to_string(),
            to_orient: '+',
            overlap: "0M".to_string(),
            tags: Vec::new(),
        });

        let (merged, mode) = merge_unambiguous_gfa(&raw);

        assert_eq!(mode, "auto_non_repeat_linear_compaction");
        assert_eq!(merged.order, vec!["a_b"]);
        assert!(merged.links.is_empty());
    }

    #[test]
    fn premerge_keeps_distinct_self_link_topologies() {
        let mut raw = test_gfa(&[("u0", "ACGT")]);
        raw.links = vec![
            Link {
                from_name: "u0".to_string(),
                from_orient: '+',
                to_name: "u0".to_string(),
                to_orient: '-',
                overlap: "0M".to_string(),
                tags: Vec::new(),
            },
            Link {
                from_name: "u0".to_string(),
                from_orient: '-',
                to_name: "u0".to_string(),
                to_orient: '+',
                overlap: "0M".to_string(),
                tags: Vec::new(),
            },
        ];

        let (merged, _) = merge_unambiguous_gfa(&raw);

        assert_eq!(merged.links.len(), 2);
    }

    #[test]
    fn premerge_keeps_conflicting_reciprocal_self_link_evidence() {
        let mut raw = test_gfa(&[("u0", "ACGT")]);
        raw.links = vec![
            Link {
                from_name: "u0".to_string(),
                from_orient: '+',
                to_name: "u0".to_string(),
                to_orient: '+',
                overlap: "3M".to_string(),
                tags: vec!["EC:i:499".to_string()],
            },
            Link {
                from_name: "u0".to_string(),
                from_orient: '-',
                to_name: "u0".to_string(),
                to_orient: '-',
                overlap: "3M".to_string(),
                tags: vec!["EC:i:500".to_string()],
            },
            Link {
                from_name: "u0".to_string(),
                from_orient: '-',
                to_name: "u0".to_string(),
                to_orient: '-',
                overlap: "2M".to_string(),
                tags: vec!["EC:i:499".to_string()],
            },
        ];

        let (merged, _) = merge_unambiguous_gfa(&raw);

        assert_eq!(merged.links.len(), 3);
    }

    fn output_guard_options(root: &Path, component: &str) -> RebuildOptions {
        fs::create_dir_all(root).unwrap();
        let root = fs::canonicalize(root).unwrap();
        let gfa = root.join("source.gfa");
        let fasta = root.join("source.fasta");
        fs::write(&gfa, "S\tnode\tACGT\n").unwrap();
        fs::write(&fasta, ">sequence\nACGT\n").unwrap();
        RebuildOptions::from_args(&[
            "--edited-gfa".into(), gfa.display().to_string(),
            "--polished-fasta".into(), fasta.display().to_string(),
            "--subgraph".into(), component.into(),
            "--out-dir".into(), root.join("06.rebuild/mito").display().to_string(),
            "--keep-debug".into(),
        ]).unwrap()
    }

    #[test]
    fn rebuild_output_guard_allows_siblings_and_scopes_force_to_component() {
        let root = temp_work_dir("rebuild-output-siblings");
        let a = output_guard_options(&root, "subgraph_001");
        let pa = OutputPaths::new(&a.out_dir, &a.subgraph);
        pa.prepare(&a).unwrap();
        fs::write(&pa.verified_gfa, "first graph").unwrap();
        fs::write(pa.debug_dir.join("path_source.fasta"), "first debug").unwrap();
        fs::write(&pa.run_report, "first report").unwrap();
        let b = output_guard_options(&root, "subgraph_001_extra");
        let pb = OutputPaths::new(&b.out_dir, &b.subgraph);
        pb.prepare(&b).unwrap();
        fs::write(&pb.verified_gfa, "second graph").unwrap();
        fs::write(pb.debug_dir.join("path_source.fasta"), "second debug").unwrap();
        fs::write(&pb.run_report, "second report").unwrap();
        assert_eq!(fs::read_to_string(&pa.verified_gfa).unwrap(), "first graph");
        assert_eq!(fs::read_to_string(&pa.run_report).unwrap(), "first report");
        assert_eq!(fs::read_to_string(pa.debug_dir.join("path_source.fasta")).unwrap(), "first debug");
        assert!(pa.prepare(&a).is_err());
        let mut force = a.clone(); force.force = true;
        pa.prepare(&force).unwrap();
        assert!(!pa.verified_gfa.exists());
        assert!(!pa.run_report.exists());
        assert!(!pa.debug_dir.join("path_source.fasta").exists());
        assert_eq!(fs::read_to_string(&pb.verified_gfa).unwrap(), "second graph");
        assert_eq!(fs::read_to_string(&pb.run_report).unwrap(), "second report");
        assert_eq!(fs::read_to_string(pb.debug_dir.join("path_source.fasta")).unwrap(), "second debug");
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn rebuild_output_guard_rejects_each_native_log_and_work_collision() {
        let root = temp_work_dir("rebuild-output-collisions");
        let options = output_guard_options(&root, "subgraph_001");
        let paths = OutputPaths::new(&options.out_dir, &options.subgraph);
        for path in paths.owned_paths(&options.subgraph) {
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(&path, "existing evidence").unwrap();
            assert!(paths.prepare(&options).is_err(), "{}", path.display());
            assert_eq!(fs::read_to_string(&path).unwrap(), "existing evidence");
            fs::remove_file(path).unwrap();
        }
        fs::create_dir_all(&paths.debug_dir).unwrap();
        assert!(paths.prepare(&options).is_err());
        assert!(!paths.subgraph_dir.exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn rebuild_output_guard_preserves_sources_even_with_force() {
        let root = temp_work_dir("rebuild-output-sources");
        let mut options = output_guard_options(&root, "subgraph_001");
        let paths = OutputPaths::new(&options.out_dir, &options.subgraph);
        paths.prepare(&options).unwrap();
        options.force = true;
        for source in [paths.verified_gfa.clone(), paths.debug_dir.join("source.fasta"), paths.run_report.clone()] {
            fs::write(&source, "protected input").unwrap();
            options.polished_fasta = source.clone();
            assert!(paths.prepare(&options).is_err());
            assert_eq!(fs::read_to_string(&source).unwrap(), "protected input");
            fs::remove_file(&source).unwrap();
        }
        // A sibling source in the shared organelle root is outside this write scope.
        let sibling = options.out_dir.join("other/source.fasta");
        fs::create_dir_all(sibling.parent().unwrap()).unwrap();
        fs::write(&sibling, "sibling input").unwrap();
        options.polished_fasta = sibling.clone();
        paths.prepare(&options).unwrap();
        assert_eq!(fs::read_to_string(sibling).unwrap(), "sibling input");
        fs::remove_dir_all(root).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn rebuild_output_guard_rejects_symlink_outputs_and_source_aliases() {
        use std::os::unix::fs::symlink;
        let root = temp_work_dir("rebuild-output-symlinks");
        let options = output_guard_options(&root, "subgraph_001");
        let paths = OutputPaths::new(&options.out_dir, &options.subgraph);
        let target = root.join("unrelated"); fs::create_dir(&target).unwrap();
        for path in [options.out_dir.clone(), paths.logs_dir.clone(),
            paths.debug_dir.parent().unwrap().to_path_buf(), paths.subgraph_dir.clone(),
            paths.debug_dir.clone(), paths.run_report.clone()] {
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            symlink(&target, &path).unwrap();
            assert!(paths.prepare(&options).is_err(), "{}", path.display());
            let mut force = options.clone(); force.force = true;
            assert!(paths.prepare(&force).is_err(), "{}", path.display());
            fs::remove_file(path).unwrap();
        }
        symlink(root.join("absent"), &paths.run_report).unwrap();
        assert!(paths.prepare(&options).is_err());
        fs::remove_file(&paths.run_report).unwrap();
        paths.prepare(&options).unwrap();
        fs::write(&paths.verified_fasta, "source through alias").unwrap();
        let alias = root.join("alias.fasta"); symlink(&paths.verified_fasta, &alias).unwrap();
        let mut force = options.clone(); force.force = true; force.polished_fasta = alias;
        assert!(paths.prepare(&force).is_err());
        assert_eq!(fs::read_to_string(&paths.verified_fasta).unwrap(), "source through alias");
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn output_paths_use_subgraph_id() {
        let paths = OutputPaths::new(Path::new("out"), "subgraph_007");
        assert_eq!(paths.debug_dir, PathBuf::from("out/debug/subgraph_007"));
        assert_eq!(
            paths.verified_gfa,
            PathBuf::from("out/subgraph_007/rebuild_subgraph_007.gfa")
        );
        assert_eq!(
            paths.verified_fasta,
            PathBuf::from("out/subgraph_007/rebuild_subgraph_007.fasta")
        );
        assert_eq!(
            paths.verified_nodes_fasta,
            PathBuf::from("out/subgraph_007/rebuild_subgraph_007_nodes.fasta")
        );
        assert_eq!(
            paths.verified_pdf,
            PathBuf::from("out/subgraph_007/rebuild_subgraph_007.pdf")
        );
        assert_eq!(
            paths.verified_svg,
            PathBuf::from("out/subgraph_007/rebuild_subgraph_007.svg")
        );
    }

    #[test]
    fn rejects_unsafe_subgraph_id() {
        assert!(parse_subgraph_id("../subgraph_001").is_err());
        assert!(parse_subgraph_id("subgraph/001").is_err());
        assert!(parse_subgraph_id("subgraph_001").is_ok());
    }

    fn test_record(sequence: &str) -> FastaRecord {
        FastaRecord {
            header: "test".to_string(),
            id: "test".to_string(),
            sequence: sequence.to_string(),
        }
    }

    fn test_gfa(nodes: &[(&str, &str)]) -> Gfa {
        let mut segments = BTreeMap::new();
        let mut order = Vec::new();
        for (name, sequence) in nodes {
            order.push((*name).to_string());
            segments.insert(
                (*name).to_string(),
                Segment {
                    name: (*name).to_string(),
                    sequence: (*sequence).to_string(),
                    tags: Vec::new(),
                },
            );
        }
        Gfa {
            headers: Vec::new(),
            segments,
            order,
            links: Vec::new(),
            other_lines: Vec::new(),
        }
    }

    fn test_mapping(hits: Vec<PafHit>) -> NodeMapping {
        NodeMapping {
            selected: hits.first().cloned(),
            accepted_hits: hits.clone(),
            all_hits: hits,
        }
    }
}



// =====================================================================
// Inlined module: rebuild_topology.rs (checked-graph topology + recorded walk)
// =====================================================================
// Checked-graph topology and the recorded resolve walk used by rebuild.
//
// This module reads explicit evidence. It never chooses a repeat arrangement
// from a FASTA sequence or changes the copy number of an unresolved repeat.



const MAX_LINEAGE_CANDIDATES: usize = 256;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct EvidenceLink {
    pub from: String,
    pub from_orientation: char,
    pub to: String,
    pub to_orientation: char,
    pub overlap: String,
    pub tags: Vec<String>,
}

impl EvidenceLink {
    /// Only exact-length simple overlaps are supported by the sequence repair.
    pub fn overlap_bases(&self) -> Result<usize, OrgraftError> {
        let digits = self
            .overlap
            .strip_suffix('M')
            .or_else(|| self.overlap.strip_suffix('='))
            .ok_or_else(|| invalid(format!("unsupported source edge overlap {}", self.overlap)))?;
        if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
            return Err(invalid(format!(
                "unsupported source edge overlap {}",
                self.overlap
            )));
        }
        digits
            .parse::<usize>()
            .map_err(|_| invalid("source edge overlap exceeds supported length"))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct LineagePart {
    pub node: String,
    pub orientation: char,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct WalkStep {
    /// A node in the original unresolved component, shared by its repeat copies.
    pub source_node: String,
    /// The literal node name in the recorded, resolved walk.
    pub instance_node: String,
    pub orientation: char,
}

#[derive(Debug, Clone)]
pub(crate) struct ResolveContract {
    pub source_sequences: BTreeMap<String, String>,
    /// Unresolved component links, with checked self-loops restored when needed.
    pub topology_links: Vec<EvidenceLink>,
    /// Links after replaying only the repeat operations recorded by resolve.
    pub resolved_links: Vec<EvidenceLink>,
    pub restored_links: Vec<EvidenceLink>,
    pub lineage: BTreeMap<String, Vec<LineagePart>>,
    pub walk: Vec<WalkStep>,
    pub circular: bool,
    pub circular_evidence: String,
    pub reference_placement_fallback: bool,
    pub alignment_reverse_complemented: bool,
    pub alignment_rotation_step: isize,
    pub path_issues: Vec<String>,
}

impl ResolveContract {
    /// Return the unique recorded edge for an oriented path transition.
    pub fn edge_between(
        &self,
        left: &WalkStep,
        right: &WalkStep,
    ) -> Result<&EvidenceLink, OrgraftError> {
        let matching: Vec<_> = self
            .resolved_links
            .iter()
            .filter(|link| {
                physical_key(link)
                    == transition_key(
                        &left.instance_node,
                        left.orientation,
                        &right.instance_node,
                        right.orientation,
                    )
            })
            .collect();
        if matching.len() != 1 {
            return Err(invalid(format!(
                "recorded walk transition {}{} -> {}{} has {} edges",
                left.instance_node,
                left.orientation,
                right.instance_node,
                right.orientation,
                matching.len()
            )));
        }
        Ok(matching[0])
    }
}

#[derive(Clone)]
struct Graph {
    sequences: BTreeMap<String, String>,
    links: Vec<EvidenceLink>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct End(String, char);
type PhysicalKey = (End, End);

fn invalid(message: impl Into<String>) -> OrgraftError {
    OrgraftError::InvalidArgument(message.into())
}
fn flip(orientation: char) -> char {
    if orientation == '+' {
        '-'
    } else {
        '+'
    }
}
fn oriented(sequence: &str, orientation: char) -> String {
    if orientation == '+' {
        sequence.to_string()
    } else {
        reverse_complement_iupac(sequence)
    }
}
fn transition_key(
    from: &str,
    from_orientation: char,
    to: &str,
    to_orientation: char,
) -> PhysicalKey {
    let a = End(from.to_string(), from_orientation);
    let b = End(to.to_string(), flip(to_orientation));
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}
fn physical_key(link: &EvidenceLink) -> PhysicalKey {
    transition_key(
        &link.from,
        link.from_orientation,
        &link.to,
        link.to_orientation,
    )
}
fn overlaps_by_key(links: &[EvidenceLink]) -> Result<BTreeMap<PhysicalKey, usize>, OrgraftError> {
    let mut overlaps = BTreeMap::new();
    for link in links {
        let overlap = link.overlap_bases()?;
        if let Some(previous) = overlaps.insert(physical_key(link), overlap) {
            if previous != overlap {
                return Err(invalid(format!(
                    "conflicting overlaps on the same physical edge {}{} -> {}{}",
                    link.from, link.from_orientation, link.to, link.to_orientation
                )));
            }
        }
    }
    Ok(overlaps)
}
fn parse_orientation(value: &str) -> Result<char, OrgraftError> {
    match value {
        "+" => Ok('+'),
        "-" => Ok('-'),
        _ => Err(invalid("invalid GFA orientation")),
    }
}
fn parse_graph(text: &str) -> Result<Graph, OrgraftError> {
    let mut sequences = BTreeMap::new();
    let mut links = Vec::new();
    for line in text.lines() {
        let f: Vec<_> = line.split('\t').collect();
        match f.first().copied().unwrap_or("") {
            "S" => {
                if f.len() < 3 || f[2] == "*" {
                    return Err(invalid(
                        "checked evidence requires concrete GFA S sequences",
                    ));
                }
                if sequences
                    .insert(f[1].to_string(), f[2].to_ascii_uppercase())
                    .is_some()
                {
                    return Err(invalid(format!("duplicate GFA node {}", f[1])));
                }
            }
            "L" => {
                if f.len() < 6 {
                    return Err(invalid("incomplete GFA L record"));
                }
                links.push(EvidenceLink {
                    from: f[1].to_string(),
                    from_orientation: parse_orientation(f[2])?,
                    to: f[3].to_string(),
                    to_orientation: parse_orientation(f[4])?,
                    overlap: f[5].to_string(),
                    tags: f[6..].iter().map(|s| s.to_string()).collect(),
                });
            }
            _ => {}
        }
    }
    if sequences.is_empty() {
        return Err(invalid("empty GFA evidence"));
    }
    if links
        .iter()
        .any(|l| !sequences.contains_key(&l.from) || !sequences.contains_key(&l.to))
    {
        return Err(invalid("GFA link references a missing node"));
    }
    Ok(Graph { sequences, links })
}

pub(crate) fn load_contract(
    checked_gfa: &Path,
    original_component_gfa: &Path,
    resolve_details: &Path,
    subgraph: &str,
) -> Result<ResolveContract, OrgraftError> {
    contract_from_text(
        &fs::read_to_string(checked_gfa)?,
        &fs::read_to_string(original_component_gfa)?,
        &fs::read_to_string(resolve_details)?,
        subgraph,
    )
}

pub(crate) fn contract_from_text(
    checked_text: &str,
    component_text: &str,
    details_text: &str,
    subgraph: &str,
) -> Result<ResolveContract, OrgraftError> {
    let checked = parse_graph(checked_text)?;
    let component = parse_graph(component_text)?;
    let details = details_for_component(details_text, subgraph)?;
    let get = |section: &str, key: &str| -> Result<&str, OrgraftError> {
        details
            .get(&(section.to_string(), key.to_string()))
            .map(String::as_str)
            .ok_or_else(|| invalid(format!("resolve details lack {section}/{subgraph}/{key}")))
    };
    let (lineage, restored_links) = verify_compaction(&checked, &component)?;
    let mut topology_links = component.links.clone();
    topology_links.extend(restored_links.iter().cloned());
    deduplicate(&mut topology_links)?;
    let selected_circular = match get("resolution", "selected_circular")? {
        "true" => true,
        "false" => false,
        _ => return Err(invalid("invalid selected_circular value")),
    };
    let selected_order = get("resolution", "selected_order")?;
    let (resolved_links, aliases) = replay_repeat_operations(&component, selected_order)?;
    // Self-loop restoration affects a singleton; it does not alter repeat A/B selection.
    let mut resolved_links = resolved_links;
    resolved_links.extend(restored_links.iter().cloned());
    deduplicate(&mut resolved_links)?;
    let mut walk = Vec::new();
    for token in get("resolution", "oriented_nodes")?.split(',') {
        let token = token.split('@').next().unwrap_or("");
        if token.len() < 2 {
            return Err(invalid("empty recorded resolve walk"));
        }
        let (instance, orientation) = token.split_at(token.len() - 1);
        let orientation = parse_orientation(orientation)?;
        let source_node = aliases.get(instance).ok_or_else(|| {
            invalid(format!(
                "recorded walk node {instance} is absent from explicit repeat operations"
            ))
        })?;
        walk.push(WalkStep {
            source_node: source_node.clone(),
            instance_node: instance.to_string(),
            orientation,
        });
    }
    if walk.is_empty() {
        return Err(invalid("empty recorded resolve walk"));
    }
    let source_nodes: BTreeSet<_> = walk.iter().map(|s| s.source_node.clone()).collect();
    if source_nodes != component.sequences.keys().cloned().collect() {
        return Err(invalid(
            "recorded resolve walk does not cover exactly the source component nodes",
        ));
    }
    let instances: BTreeSet<_> = walk.iter().map(|s| s.instance_node.clone()).collect();
    if instances.len() != walk.len() || instances != aliases.keys().cloned().collect() {
        return Err(invalid(
            "recorded resolve walk omits or repeats a resolved node instance",
        ));
    }
    // This additional circle source is restricted to an explicit checked singleton self-loop.
    let checked_singleton_circle = lineage.len() == 1
        && lineage.values().next().is_some_and(|x| x.len() == 1)
        && topology_links
            .iter()
            .any(|l| l.from == l.to && l.from_orientation == l.to_orientation);
    // Reference-placement fallback sometimes omitted the circular flag even when
    // the recorded walk already has every edge of an explicit directed cycle.
    // This uses graph edges only, never a circular match to the target FASTA.
    let unique_transition = |left: &WalkStep, right: &WalkStep| {
        resolved_links
            .iter()
            .filter(|edge| {
                physical_key(edge)
                    == transition_key(
                        &left.instance_node,
                        left.orientation,
                        &right.instance_node,
                        right.orientation,
                    )
            })
            .count()
            == 1
    };
    let explicit_closed_walk = walk
        .windows(2)
        .all(|pair| unique_transition(&pair[0], &pair[1]))
        && unique_transition(walk.last().unwrap(), &walk[0]);
    let circular = selected_circular || checked_singleton_circle || explicit_closed_walk;
    let circular_evidence = if selected_circular {
        "recorded_resolve_selected_circular"
    } else if checked_singleton_circle {
        "explicit_checked_singleton_self_loop"
    } else if explicit_closed_walk {
        "explicit_graph_closed_recorded_walk"
    } else {
        "no_explicit_circle_evidence"
    }
    .to_string();
    let alignment_reverse_complemented = match get("alignment", "reverse_complemented")? {
        "true" => true,
        "false" => false,
        _ => return Err(invalid("invalid reverse_complemented value")),
    };
    let alignment_rotation_step = get("alignment", "rotation_step")?
        .parse::<isize>()
        .map_err(|_| invalid("invalid recorded alignment rotation"))?;
    let mut contract = ResolveContract {
        source_sequences: component.sequences,
        topology_links,
        resolved_links,
        restored_links,
        lineage,
        walk,
        circular,
        circular_evidence,
        reference_placement_fallback: get("resolution", "score_method")?
            == "reference-placement-fallback",
        alignment_reverse_complemented,
        alignment_rotation_step,
        path_issues: Vec::new(),
    };
    for pair in contract.walk.windows(2) {
        if let Err(error) = contract.edge_between(&pair[0], &pair[1]) {
            contract.path_issues.push(error.to_string());
        }
    }
    if contract.circular {
        if let Err(error) = contract.edge_between(contract.walk.last().unwrap(), &contract.walk[0])
        {
            contract.path_issues.push(error.to_string());
        }
    }
    Ok(contract)
}

fn details_for_component(
    text: &str,
    subgraph: &str,
) -> Result<BTreeMap<(String, String), String>, OrgraftError> {
    let mut rows = BTreeMap::new();
    let mut lines = text.lines();
    if lines.next().unwrap_or("") != "section\tid\tkey\tvalue\tnotes" {
        return Err(invalid(
            "resolve details header does not match the recorded TSV contract",
        ));
    }
    for line in lines {
        let f: Vec<_> = line.splitn(5, '\t').collect();
        if f.len() < 4 || f[1] != subgraph || !matches!(f[0], "resolution" | "alignment") {
            continue;
        }
        let key = (f[0].to_string(), f[2].to_string());
        if let Some(previous) = rows.insert(key, f[3].to_string()) {
            if previous != f[3] {
                return Err(invalid("conflicting resolve detail rows"));
            }
        }
    }
    Ok(rows)
}

fn tokenize_name(
    name: &str,
    checked: &Graph,
    prefix: &mut Vec<String>,
    output: &mut Vec<Vec<String>>,
) -> Result<(), OrgraftError> {
    if output.len() >= MAX_LINEAGE_CANDIDATES {
        return Err(invalid("too many possible linear-compaction name mappings"));
    }
    for node in checked.sequences.keys() {
        if name == node {
            let mut path = prefix.clone();
            path.push(node.clone());
            output.push(path);
        } else if let Some(rest) = name.strip_prefix(&format!("{node}_")) {
            prefix.push(node.clone());
            tokenize_name(rest, checked, prefix, output)?;
            prefix.pop();
        }
    }
    Ok(())
}
fn signed_paths(
    nodes: &[String],
    index: usize,
    offset: usize,
    target: &str,
    checked: &Graph,
    path: &mut Vec<LineagePart>,
    output: &mut Vec<Vec<LineagePart>>,
) -> Result<(), OrgraftError> {
    if output.len() >= MAX_LINEAGE_CANDIDATES {
        return Err(invalid("too many possible signed compaction paths"));
    }
    if index == nodes.len() {
        if offset == target.len() {
            output.push(path.clone());
        }
        return Ok(());
    }
    let node = &nodes[index];
    for orientation in ['+', '-'] {
        let sequence = oriented(&checked.sequences[node], orientation);
        if !target
            .get(offset..)
            .is_some_and(|rest| rest.starts_with(&sequence))
        {
            continue;
        }
        if let Some(previous) = path.last() {
            let key = transition_key(&previous.node, previous.orientation, node, orientation);
            if !checked.links.iter().any(|l| physical_key(l) == key) {
                continue;
            }
        }
        path.push(LineagePart {
            node: node.clone(),
            orientation,
        });
        signed_paths(
            nodes,
            index + 1,
            offset + sequence.len(),
            target,
            checked,
            path,
            output,
        )?;
        path.pop();
    }
    Ok(())
}

/// Reuse the checked-graph compaction proof for a derived fixed walk.
/// The mapping is unique by source sequence, direction and every physical edge.
pub(crate) fn verified_compaction_lineage(
    checked: &Path,
    component: &Path,
) -> Result<BTreeMap<String, Vec<LineagePart>>, OrgraftError> {
    let checked = parse_graph(&fs::read_to_string(checked)?)?;
    let component = parse_graph(&fs::read_to_string(component)?)?;
    let (lineage, restored) = verify_compaction(&checked, &component)?;
    let nodes = lineage
        .values()
        .flatten()
        .map(|part| part.node.clone())
        .collect::<Vec<_>>();
    if nodes.len() != checked.sequences.len()
        || nodes.iter().cloned().collect::<BTreeSet<_>>()
            != checked.sequences.keys().cloned().collect()
    {
        return Err(invalid(
            "derived compaction must cover every checked node exactly once",
        ));
    }
    if !restored.is_empty() {
        return Err(invalid(
            "derived compaction must retain every physical edge explicitly",
        ));
    }
    Ok(lineage)
}

fn verify_compaction(
    checked: &Graph,
    component: &Graph,
) -> Result<(BTreeMap<String, Vec<LineagePart>>, Vec<EvidenceLink>), OrgraftError> {
    let mut candidates = Vec::new();
    for (name, sequence) in &component.sequences {
        let possible = if let Some(original) = checked.sequences.get(name) {
            if original != sequence {
                return Err(invalid(format!(
                    "source component node {name} differs from checked evidence"
                )));
            }
            vec![vec![LineagePart {
                node: name.clone(),
                orientation: '+',
            }]]
        } else {
            let mut tokenized = Vec::new();
            tokenize_name(name, checked, &mut Vec::new(), &mut tokenized)?;
            let mut possible = Vec::new();
            for nodes in tokenized {
                signed_paths(
                    &nodes,
                    0,
                    0,
                    sequence,
                    checked,
                    &mut Vec::new(),
                    &mut possible,
                )?;
            }
            possible
        };
        if possible.is_empty() {
            return Err(invalid(format!(
                "no sequence-and-link-supported checked lineage for {name}"
            )));
        }
        candidates.push((name.clone(), possible));
    }
    let mut accepted = Vec::new();
    let mut last_rejection = None;
    try_lineages(
        &candidates,
        0,
        &mut BTreeMap::new(),
        checked,
        component,
        &mut accepted,
        &mut 0,
        &mut last_rejection,
    )?;
    if accepted.len() != 1 {
        return Err(invalid(format!(
            "checked-to-component lineage has {} valid interpretations{}",
            accepted.len(),
            last_rejection
                .map(|reason| format!(": {reason}"))
                .unwrap_or_default()
        )));
    }
    Ok(accepted.remove(0))
}
fn try_lineages(
    candidates: &[(String, Vec<Vec<LineagePart>>)],
    index: usize,
    current: &mut BTreeMap<String, Vec<LineagePart>>,
    checked: &Graph,
    component: &Graph,
    accepted: &mut Vec<(BTreeMap<String, Vec<LineagePart>>, Vec<EvidenceLink>)>,
    attempts: &mut usize,
    last_rejection: &mut Option<String>,
) -> Result<(), OrgraftError> {
    if accepted.len() > 1 {
        return Ok(());
    }
    if index == candidates.len() {
        *attempts += 1;
        if *attempts > MAX_LINEAGE_CANDIDATES {
            return Err(invalid("checked compaction mapping search limit reached"));
        }
        match projected_links(checked, component, current) {
            Ok(restored) => accepted.push((current.clone(), restored)),
            Err(error) => *last_rejection = Some(error.to_string()),
        }
        return Ok(());
    }
    let (name, choices) = &candidates[index];
    for choice in choices {
        current.insert(name.clone(), choice.clone());
        try_lineages(
            candidates,
            index + 1,
            current,
            checked,
            component,
            accepted,
            attempts,
            last_rejection,
        )?;
        current.remove(name);
    }
    Ok(())
}
fn projected_links(
    checked: &Graph,
    component: &Graph,
    lineage: &BTreeMap<String, Vec<LineagePart>>,
) -> Result<Vec<EvidenceLink>, OrgraftError> {
    let mut membership = BTreeMap::new();
    for (merged, parts) in lineage {
        for (index, part) in parts.iter().enumerate() {
            if membership
                .insert(
                    part.node.clone(),
                    (merged.clone(), index, parts.len(), part.orientation),
                )
                .is_some()
            {
                return Err(invalid("checked node used by more than one compacted node"));
            }
        }
    }
    let endpoint = |node: &str, side: char| -> Option<(String, char)> {
        let (merged, index, count, orientation) = membership.get(node)?;
        let side = if *orientation == '+' {
            side
        } else {
            flip(side)
        };
        if *index == 0 && side == '-' {
            Some((merged.clone(), '-'))
        } else if *index + 1 == *count && side == '+' {
            Some((merged.clone(), '+'))
        } else {
            None
        }
    };
    let mut projected = Vec::new();
    let mut internal = BTreeSet::new();
    for link in &checked.links {
        let left = membership.get(&link.from);
        let right = membership.get(&link.to);
        match (left, right) {
            (None, None) => continue,
            (None, Some(_)) | (Some(_), None) => {
                return Err(invalid("source component excludes a linked checked node"))
            }
            _ => {}
        }
        let (lm, li, lc, lo) = left.unwrap();
        let (rm, ri, _, ro) = right.unwrap();
        if lm == rm && *lc > 1 {
            if li.abs_diff(*ri) != 1 {
                return Err(invalid(
                    "unexpected internal link in a proposed linear compaction",
                ));
            }
            let expected = if li < ri {
                transition_key(&link.from, *lo, &link.to, *ro)
            } else {
                transition_key(&link.to, *ro, &link.from, *lo)
            };
            if physical_key(link) != expected {
                return Err(invalid(
                    "internal compaction edge has incompatible orientations",
                ));
            }
            internal.insert(if link.from < link.to {
                (link.from.clone(), link.to.clone())
            } else {
                (link.to.clone(), link.from.clone())
            });
            continue;
        }
        let (from, from_orientation) = endpoint(&link.from, link.from_orientation)
            .ok_or_else(|| invalid("external link enters an internal compacted endpoint"))?;
        let (to, to_side) = endpoint(&link.to, flip(link.to_orientation))
            .ok_or_else(|| invalid("external link enters an internal compacted endpoint"))?;
        projected.push(EvidenceLink {
            from,
            from_orientation,
            to,
            to_orientation: flip(to_side),
            overlap: link.overlap.clone(),
            tags: link.tags.clone(),
        });
    }
    for parts in lineage.values() {
        for pair in parts.windows(2) {
            let mut names = [pair[0].node.clone(), pair[1].node.clone()];
            names.sort();
            if !internal.contains(&(names[0].clone(), names[1].clone())) {
                return Err(invalid("missing internal linear-compaction edge"));
            }
        }
    }
    let expected_overlaps = overlaps_by_key(&projected)?;
    let existing_overlaps = overlaps_by_key(&component.links)?;
    for (key, overlap) in &existing_overlaps {
        match expected_overlaps.get(key) {
            None => return Err(invalid("component has links absent from checked topology")),
            Some(expected) if expected != overlap => {
                return Err(invalid(
                    "component edge overlap differs from checked evidence",
                ))
            }
            Some(_) => {}
        }
    }
    let expected: BTreeSet<_> = expected_overlaps.keys().cloned().collect();
    let existing: BTreeSet<_> = existing_overlaps.keys().cloned().collect();
    let missing: BTreeSet<_> = expected.difference(&existing).cloned().collect();
    let mut restored = Vec::new();
    let mut seen = BTreeSet::new();
    for link in projected {
        let key = physical_key(&link);
        if !missing.contains(&key) || !seen.insert(key) {
            continue;
        }
        // A single checked node is a complete component here: the membership
        // check above already rejected links to every other checked node.
        if membership.len() != 1
            || link.from != link.to
            || link.from_orientation != link.to_orientation
        {
            return Err(invalid(
                "checked topology differs beyond a missing singleton self-loop",
            ));
        }
        restored.push(link);
    }
    Ok(restored)
}

fn replay_repeat_operations(
    component: &Graph,
    operations: &str,
) -> Result<(Vec<EvidenceLink>, BTreeMap<String, String>), OrgraftError> {
    let mut aliases: BTreeMap<_, _> = component
        .sequences
        .keys()
        .map(|n| (n.clone(), n.clone()))
        .collect();
    let mut links = component.links.clone();
    // Reverse-written records are one physical edge, even without repeat edits.
    // Keep the first record and its tags; conflicting overlaps still fail.
    deduplicate(&mut links)?;
    if operations == "." || operations.is_empty() {
        return Ok((links, aliases));
    }
    for operation in operations.split(',') {
        let (node, action) = operation
            .rsplit_once(':')
            .ok_or_else(|| invalid("malformed recorded repeat operation"))?;
        let (strategy, duplicate) = action
            .split_once("->")
            .ok_or_else(|| invalid("malformed recorded repeat copy"))?;
        if !matches!(strategy, "A" | "B") || aliases.contains_key(duplicate) {
            return Err(invalid("invalid recorded repeat operation"));
        }
        let source = aliases
            .get(node)
            .cloned()
            .ok_or_else(|| invalid("recorded repeat source does not exist"))?;
        aliases.insert(duplicate.to_string(), source);
        let indices = |side: char| -> Vec<usize> {
            links
                .iter()
                .enumerate()
                .filter_map(|(i, l)| {
                    let endpoint = if l.from == node && l.to == node {
                        None
                    } else if l.from == node {
                        Some(l.from_orientation)
                    } else if l.to == node {
                        Some(flip(l.to_orientation))
                    } else {
                        None
                    };
                    (endpoint == Some(side)).then_some(i)
                })
                .collect()
        };
        let minus = indices('-');
        let plus = indices('+');
        if minus.len() != 2 || plus.len() != 2 {
            return Err(invalid(
                "recorded repeat operation no longer has two links on each end",
            ));
        }
        let incident: Vec<_> = links
            .iter()
            .enumerate()
            .filter_map(|(i, l)| ((l.from == node) || (l.to == node)).then_some(i))
            .collect();
        let mut copied = BTreeMap::new();
        for &index in &incident {
            let mut link = links[index].clone();
            if link.from == node {
                link.from = duplicate.to_string();
            }
            if link.to == node {
                link.to = duplicate.to_string();
            }
            copied.insert(index, links.len());
            links.push(link);
        }
        let plus_index = if strategy == "A" { 0 } else { 1 };
        let keep = BTreeSet::from([
            minus[0],
            plus[plus_index],
            copied[&minus[1]],
            copied[&plus[1 - plus_index]],
        ]);
        let affected: BTreeSet<_> = incident
            .into_iter()
            .chain(copied.values().copied())
            .collect();
        links = links
            .into_iter()
            .enumerate()
            .filter_map(|(i, l)| (!affected.contains(&i) || keep.contains(&i)).then_some(l))
            .collect();
        deduplicate(&mut links)?;
    }
    Ok((links, aliases))
}
fn deduplicate(links: &mut Vec<EvidenceLink>) -> Result<(), OrgraftError> {
    overlaps_by_key(links)?;
    let mut seen = BTreeSet::new();
    links.retain(|link| seen.insert(physical_key(link)));
    Ok(())
}

#[cfg(test)]
mod tests_topology {
    use super::*;
    fn details(walk: &str, circular: bool, operations: &str) -> String {
        format!("section\tid\tkey\tvalue\tnotes\nresolution\tsubgraph_001\tselected_order\t{operations}\t\nresolution\tsubgraph_001\toriented_nodes\t{walk}\t\nresolution\tsubgraph_001\tselected_circular\t{circular}\t\nresolution\tsubgraph_001\tscore_method\treference-placement-fallback\t\nalignment\tsubgraph_001\treverse_complemented\tfalse\t\nalignment\tsubgraph_001\trotation_step\t0\t\n")
    }
    #[test]
    fn restores_only_explicit_checked_singleton_self_loop() {
        let checked = "S\ta\tACGT\nL\ta\t+\ta\t+\t0M\n";
        let c = contract_from_text(
            checked,
            "S\ta\tACGT\n",
            &details("a+", false, "."),
            "subgraph_001",
        )
        .unwrap();
        assert!(c.circular);
        assert_eq!(c.restored_links.len(), 1);
        assert!(c.path_issues.is_empty());
        let linear = contract_from_text(
            "S\ta\tACGT\n",
            "S\ta\tACGT\n",
            &details("a+", false, "."),
            "subgraph_001",
        )
        .unwrap();
        assert!(!linear.circular);
    }
    #[test]
    fn rejects_component_that_omits_a_checked_neighbor() {
        let checked = "S\ta\tAAA\nS\tb\tCCC\nL\ta\t+\tb\t+\t0M\n";
        assert!(contract_from_text(
            checked,
            "S\ta\tAAA\n",
            &details("a+", false, "."),
            "subgraph_001"
        )
        .is_err());
    }
    #[test]
    fn explicit_two_node_cycle_survives_a_missing_fallback_circle_flag() {
        let graph = "S\ta\tACG\nS\tb\tTTA\nL\ta\t+\tb\t+\t0M\nL\ta\t-\tb\t-\t0M\n";
        let c = contract_from_text(graph, graph, &details("b+,a+", false, "."), "subgraph_001")
            .unwrap();
        assert!(c.circular);
        assert_eq!(c.circular_evidence, "explicit_graph_closed_recorded_walk");
        assert!(c.restored_links.is_empty());
        assert!(c.path_issues.is_empty());
    }
    #[test]
    fn a_closing_edge_does_not_excuse_an_invalid_internal_walk() {
        let graph = "S\ta\tACG\nS\tb\tTTA\nS\tc\tGGC\nL\ta\t+\tb\t+\t0M\nL\tb\t+\tc\t+\t0M\nL\tc\t+\ta\t+\t0M\nL\tb\t+\ta\t+\t0M\n";
        let c = contract_from_text(
            graph,
            graph,
            &details("a+,c+,b+", false, "."),
            "subgraph_001",
        )
        .unwrap();
        assert!(!c.circular);
        assert!(!c.path_issues.is_empty());
    }
    #[test]
    fn keeps_invalid_reference_placement_explicit() {
        let graph = "S\ta\tAAA\nS\tb\tCCC\nL\ta\t+\tb\t-\t0M\n";
        let c = contract_from_text(graph, graph, &details("a+,b+", false, "."), "subgraph_001")
            .unwrap();
        assert_eq!(c.path_issues.len(), 1);
    }
    #[test]
    fn matches_edges_by_direction_even_with_different_overlaps() {
        let graph = "S\ta\tAAAAC\nS\tb\tACCCC\nL\ta\t+\tb\t+\t2M\nL\ta\t-\tb\t-\t3M\n";
        let c = contract_from_text(graph, graph, &details("b+,a+", false, "."), "subgraph_001")
            .unwrap();
        assert_eq!(
            c.edge_between(&c.walk[0], &c.walk[1]).unwrap().overlap,
            "3M"
        );
    }

    #[test]
    fn rejects_unrecorded_component_overlap_change() {
        let checked = "S\ta\tAAC\nS\tb\tACC\nL\ta\t+\tb\t+\t2M\n";
        let component = "S\ta\tAAC\nS\tb\tACC\nL\ta\t+\tb\t+\t0M\n";
        let error = contract_from_text(
            checked,
            component,
            &details("a+,b+", false, "."),
            "subgraph_001",
        )
        .unwrap_err();
        assert!(error
            .to_string()
            .contains("component edge overlap differs from checked evidence"));
    }

    #[test]
    fn conflicting_parallel_overlap_is_rejected_before_repeat_deduplication() {
        let first = EvidenceLink {
            from: "a".into(),
            from_orientation: '+',
            to: "r".into(),
            to_orientation: '+',
            overlap: "0M".into(),
            tags: Vec::new(),
        };
        let mut second = first.clone();
        second.overlap = "1M".into();
        let mut links = vec![first, second];
        let error = deduplicate(&mut links).unwrap_err();
        assert!(error.to_string().contains("conflicting overlaps"));
        assert_eq!(links.len(), 2);
    }

    #[test]
    fn reverse_duplicate_with_equal_overlap_keeps_the_same_evidence() {
        let mut graph =
            parse_graph("S\ta\tAAC\nS\tb\tACC\nL\ta\t+\tb\t+\t2M\nL\tb\t-\ta\t-\t2=\n").unwrap();
        deduplicate(&mut graph.links).unwrap();
        assert_eq!(graph.links.len(), 1);
        assert_eq!(graph.links[0].overlap, "2M");
    }

    #[test]
    fn contract_deduplicates_physical_self_edges_with_or_without_operations() {
        let first = "L\ta\t+\ta\t+\t0M\tRC:i:26\tSC:Z:first\n";
        let reverse = "L\ta\t-\ta\t-\t0=\tRC:i:14\tSC:Z:second\n";
        let checked = format!("S\ta\tACGT\n{first}{reverse}");
        for operations in [".", ""] {
            for component in [checked.as_str(), "S\ta\tACGT\n"] {
                let contract = contract_from_text(
                    &checked,
                    component,
                    &details("a+", true, operations),
                    "subgraph_001",
                )
                .unwrap();
                assert_eq!(contract.topology_links.len(), 1);
                assert_eq!(contract.resolved_links.len(), 1);
                let link = contract
                    .edge_between(&contract.walk[0], &contract.walk[0])
                    .unwrap();
                assert_eq!(link.from_orientation, '+');
                assert_eq!(link.to_orientation, '+');
                assert_eq!(link.overlap, "0M");
                assert_eq!(link.tags, vec!["RC:i:26", "SC:Z:first"]);
                assert_eq!(contract.topology_links[0].tags, link.tags);
            }
        }
    }

    #[test]
    fn contract_physical_edge_dedup_rejects_conflicts_and_keeps_distinct_edges() {
        let conflicting = "S\ta\tACGT\nL\ta\t+\ta\t+\t0M\nL\ta\t-\ta\t-\t1M\n";
        for component in [conflicting, "S\ta\tACGT\n"] {
            let error = contract_from_text(
                conflicting,
                component,
                &details("a+", true, "."),
                "subgraph_001",
            )
            .err()
            .unwrap();
            assert!(error.to_string().contains("conflicting overlaps"));
        }
        let graph = "S\ta\tACGT\nS\tb\tTGCA\nL\ta\t+\tb\t+\t0M\tRC:i:7\nL\tb\t-\ta\t-\t0M\tRC:i:9\nL\ta\t+\tb\t-\t0M\tRC:i:11\n";
        let contract =
            contract_from_text(graph, graph, &details("a+,b+", false, "."), "subgraph_001").unwrap();
        assert_eq!(contract.topology_links.len(), 2);
        assert_eq!(contract.resolved_links.len(), 2);
        assert_eq!(contract.topology_links[0].tags, vec!["RC:i:7"]);
        assert_eq!(contract.topology_links[1].tags, vec!["RC:i:11"]);
    }

    #[test]
    fn unrelated_checked_component_overlap_is_outside_the_requested_contract() {
        let component = "S\ta\tAAC\nS\tb\tACC\nL\ta\t+\tb\t+\t2M\n";
        let checked = format!("{component}S\tz\tTTT\nL\tz\t+\tz\t+\t*\n");
        let contract = contract_from_text(
            &checked,
            component,
            &details("a+,b+", false, "."),
            "subgraph_001",
        )
        .unwrap();
        assert!(contract.path_issues.is_empty());
    }
}



// =====================================================================
// Inlined module: rebuild_path_alignment.rs (fixed-spell -> FASTA alignment)
// =====================================================================
// Complete alignments from a fixed graph spell to an authoritative FASTA.
//
// Only the recorded overall direction and, for a circle, a coordinate rotation
// may change. A partial minimap2 hit is an anchor, never a complete result.


const MAX_TERMINAL_QUERY_BP: usize = 2_048;
const MAX_TERMINAL_TARGET_BP: usize = 4_096;
const MAX_DP_CELLS: usize = 4_500_000;
const MAX_CANDIDATES: usize = 128;
const MAX_TRACEBACK_STATES: usize = 100_000;
// A complete global alignment has one canonical traceback per supported frame.
// Its memory and work limits describe the implementation, not graph validity.
const MAX_GLOBAL_WAVEFRONT_CELLS: usize = 4_500_000;
const MAX_GLOBAL_MATCH_WORK: usize = 50_000_000;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct AlignmentCandidate {
    pub target_in_path_frame: String,
    /// Offset in target after applying the explicitly recorded overall strand.
    pub frame_offset: usize,
    /// Source is query: I consumes source, D consumes target.
    pub cigar: String,
    pub method: String,
    /// An evidenced repeated target prefix after one full circle, never clipped query.
    pub terminal_target_overlap: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Op {
    count: usize,
    code: u8,
}

#[derive(Debug)]
struct Anchor {
    query_start: usize,
    query_end: usize,
    target_start: usize,
    target_end: usize,
    cigar: Vec<Op>,
    mapq: usize,
    primary: bool,
}

fn reverse_complement_iupac(sequence: &str) -> String {
    sequence
        .chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'R' => 'Y',
            'Y' => 'R',
            'M' => 'K',
            'K' => 'M',
            'B' => 'V',
            'V' => 'B',
            'D' => 'H',
            'H' => 'D',
            other => other,
        })
        .collect()
}
fn dna_iupac_string(sequence: &str) -> Result<String, String> {
    let value = sequence.to_ascii_uppercase();
    if value.is_empty() || !value.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b)) {
        return Err("alignment_input_is_empty_or_not_iupac_dna".to_string());
    }
    Ok(value)
}
fn rotate(sequence: &str, offset: usize) -> String {
    let offset = offset % sequence.len();
    format!("{}{}", &sequence[offset..], &sequence[..offset])
}
fn periodic_slice(sequence: &[u8], start: usize, length: usize) -> Vec<u8> {
    (0..length)
        .map(|i| sequence[(start + i) % sequence.len()])
        .collect()
}
fn push_op(ops: &mut Vec<Op>, count: usize, code: u8) {
    if count == 0 {
        return;
    }
    if let Some(last) = ops.last_mut() {
        if last.code == code {
            last.count += count;
            return;
        }
    }
    ops.push(Op { count, code });
}
fn append_ops(target: &mut Vec<Op>, more: &[Op]) {
    for op in more {
        push_op(target, op.count, op.code);
    }
}
fn encode(ops: &[Op]) -> String {
    ops.iter()
        .map(|op| format!("{}{}", op.count, op.code as char))
        .collect()
}
fn parse_cigar(cigar: &str) -> Result<Vec<Op>, String> {
    let mut ops = Vec::new();
    let mut count = 0usize;
    let mut digits = false;
    for code in cigar.bytes() {
        if code.is_ascii_digit() {
            count = count
                .checked_mul(10)
                .and_then(|n| n.checked_add((code - b'0') as usize))
                .ok_or("alignment_cigar_length_overflow")?;
            digits = true;
        } else {
            if !digits || count == 0 || !b"M=XID".contains(&code) {
                return Err("alignment_cigar_has_unsupported_operation".to_string());
            }
            push_op(&mut ops, count, code);
            count = 0;
            digits = false;
        }
    }
    if digits || ops.is_empty() {
        return Err("invalid_alignment_cigar".to_string());
    }
    Ok(ops)
}
fn consumed(ops: &[Op]) -> Result<(usize, usize), String> {
    let mut q = 0usize;
    let mut t = 0usize;
    for op in ops {
        if op.code != b'D' {
            q = q.checked_add(op.count).ok_or("alignment_length_overflow")?;
        }
        if op.code != b'I' {
            t = t.checked_add(op.count).ok_or("alignment_length_overflow")?;
        }
    }
    Ok((q, t))
}
fn normalized_cigar(source: &[u8], target: &[u8], ops: &[Op]) -> Result<String, String> {
    if consumed(ops)? != (source.len(), target.len()) {
        return Err("alignment_does_not_consume_both_complete_sequences".to_string());
    }
    let mut q = 0;
    let mut t = 0;
    let mut normalized = Vec::new();
    for op in ops {
        match op.code {
            b'M' | b'=' | b'X' => {
                for _ in 0..op.count {
                    let equal = source[q] == target[t];
                    if op.code == b'=' && !equal {
                        return Err("alignment_equality_conflicts_with_sequence".to_string());
                    }
                    if op.code == b'X' && equal {
                        return Err(
                            "alignment_mismatch_operator_conflicts_with_sequence".to_string()
                        );
                    }
                    push_op(&mut normalized, 1, if equal { b'=' } else { b'X' });
                    q += 1;
                    t += 1;
                }
            }
            b'I' => {
                push_op(&mut normalized, op.count, b'I');
                q += op.count;
            }
            b'D' => {
                push_op(&mut normalized, op.count, b'D');
                t += op.count;
            }
            _ => return Err("unsupported_global_alignment_operation".to_string()),
        }
    }
    Ok(encode(&normalized))
}

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|error| {
            format!(
                "cannot_create_alignment_evidence {}: {error}",
                path.display()
            )
        })?;
    file.write_all(bytes)
        .map_err(|error| format!("cannot_write_alignment_evidence: {error}"))
}

pub(crate) fn align_path(
    source: &str,
    target: &str,
    circular: bool,
    target_reverse: bool,
    minimap2: &Path,
    threads: usize,
    work_dir: &Path,
) -> Result<Vec<AlignmentCandidate>, String> {
    let source = dna_iupac_string(source)?;
    let target = dna_iupac_string(target)?;
    let target = if target_reverse {
        reverse_complement_iupac(&target)
    } else {
        target
    };
    if source == target {
        return Ok(vec![AlignmentCandidate {
            cigar: format!("{}=", source.len()),
            target_in_path_frame: target,
            frame_offset: 0,
            method: "exact_same_frame".to_string(),
            terminal_target_overlap: 0,
        }]);
    }
    if circular && source.len() == target.len() {
        if let Some(offset) = format!("{target}{target}").find(&source) {
            if offset < target.len() {
                return Ok(vec![AlignmentCandidate {
                    target_in_path_frame: source.clone(),
                    frame_offset: offset,
                    cigar: format!("{}=", source.len()),
                    method: "exact_circular_frame".to_string(),
                    terminal_target_overlap: 0,
                }]);
            }
        }
    }
    if threads == 0 {
        return Err("alignment_threads_must_be_positive".to_string());
    }
    fs::create_dir_all(work_dir)
        .map_err(|error| format!("cannot_create_alignment_work_directory: {error}"))?;
    let source_path = work_dir.join("path_source.fasta");
    let target_path = work_dir.join("path_target.fasta");
    write_new(&source_path, format!(">path_source\n{source}\n").as_bytes())?;
    let target_for_alignment = if circular {
        format!("{target}{target}")
    } else {
        target.clone()
    };
    write_new(
        &target_path,
        format!(">path_target\n{target_for_alignment}\n").as_bytes(),
    )?;
    let output = Command::new(minimap2)
        .args(["-x", "asm5", "-c", "--eqx", "--for-only", "-N", "100", "-t"])
        .arg(threads.to_string())
        .arg(&target_path)
        .arg(&source_path)
        .output()
        .map_err(|error| format!("cannot_run_path_minimap2: {error}"))?;
    write_new(&work_dir.join("path_minimap2.paf"), &output.stdout)?;
    write_new(&work_dir.join("path_minimap2.stderr.log"), &output.stderr)?;
    if !output.status.success() {
        return Err(format!(
            "path_minimap2_failed_exit_{:?}",
            output.status.code()
        ));
    }
    let paf = String::from_utf8(output.stdout).map_err(|_| "path_minimap2_output_is_not_utf8")?;
    match candidates_from_paf(&source, &target, circular, &paf) {
        Ok(candidates) => Ok(candidates),
        Err(reason) if can_try_complete_global_alignment(&reason) => {
            let result = global_candidates_from_paf(&source, &target, circular, &paf);
            let mut evidence = format!("previous_bounded_alignment\t{reason}\n");
            match &result {
                Ok(candidates) => {
                    for candidate in candidates {
                        evidence.push_str(&format!(
                            "complete_global_candidate\t{}\t{}\t{}\n",
                            candidate.frame_offset, candidate.cigar, candidate.method
                        ));
                    }
                }
                Err(error) => evidence.push_str(&format!("global_alignment_unresolved\t{error}\n")),
            }
            write_new(
                &work_dir.join("path_global_alignment.tsv"),
                evidence.as_bytes(),
            )?;
            result.map_err(|error| format!("{error}; previous_bounded_alignment={reason}"))
        }
        Err(reason) => Err(reason),
    }
}

fn can_try_complete_global_alignment(reason: &str) -> bool {
    // Never replace malformed CIGAR or strand evidence with a convenient alignment.
    if reason.contains("conflicts_with_sequence")
        || reason.contains("operator_conflicts_with_sequence")
        || reason.contains("paf_")
        || reason.contains("not_a_forward")
    {
        return false;
    }
    [
        "long_source_ends",
        "long_target_ends",
        "more_than_one_target_copy",
        "resource_exhausted",
        "terminal_extension_size_overflow",
    ]
    .iter()
    .any(|kind| reason.contains(kind))
}

struct WavefrontLayer {
    /// Furthest query coordinate after extending exact matches on each diagonal.
    query: Vec<isize>,
    transition: Vec<u8>,
}

impl WavefrontLayer {
    fn at(&self, score: usize, diagonal: isize) -> Option<isize> {
        if diagonal.unsigned_abs() > score {
            return None;
        }
        let value = self.query[(diagonal + score as isize) as usize];
        (value >= 0).then_some(value)
    }
}

struct GlobalAlignment {
    cigar: String,
    edit_distance: usize,
    explored_cells: usize,
}

/// Exact unit-cost global edit alignment, with a deterministic optimal traceback.
/// It does not enumerate equivalent gap placements. Shared-copy boundary choices
/// remain the subsequent projection module's job; every returned CIGAR is full.
fn complete_global_cigar(source: &[u8], target: &[u8]) -> Result<GlobalAlignment, String> {
    complete_global_cigar_with_limits(
        source,
        target,
        MAX_GLOBAL_WAVEFRONT_CELLS,
        MAX_GLOBAL_MATCH_WORK,
    )
}

fn complete_global_cigar_with_limits(
    source: &[u8],
    target: &[u8],
    max_cells: usize,
    max_work: usize,
) -> Result<GlobalAlignment, String> {
    let m = isize::try_from(source.len()).map_err(|_| "global_sequence_length_overflow")?;
    let n = isize::try_from(target.len()).map_err(|_| "global_sequence_length_overflow")?;
    let end_diagonal = m - n;
    let mut work = 0usize;
    let mut extend = |mut q: isize, k: isize| -> Result<isize, String> {
        let mut t = q - k;
        while q < m && t < n {
            work += 1;
            if work > max_work {
                return Err("global_wavefront_work_limit".to_string());
            }
            if source[q as usize] != target[t as usize] {
                break;
            }
            q += 1;
            t += 1;
        }
        Ok(q)
    };
    let first = extend(0, 0)?;
    let mut layers = vec![WavefrontLayer {
        query: vec![first],
        transition: vec![b'S'],
    }];
    let mut score = 0usize;
    let mut explored_cells = 1usize;
    loop {
        if layers[score].at(score, end_diagonal) == Some(m) {
            break;
        }
        score += 1;
        let width = score
            .checked_mul(2)
            .and_then(|v| v.checked_add(1))
            .ok_or("global_wavefront_size_overflow")?;
        explored_cells = explored_cells
            .checked_add(width)
            .ok_or("global_wavefront_size_overflow")?;
        if explored_cells > max_cells {
            return Err("global_wavefront_memory_limit".to_string());
        }
        let previous = &layers[score - 1];
        let mut layer = WavefrontLayer {
            query: vec![-1; width],
            transition: vec![0; width],
        };
        for k in -(score as isize)..=score as isize {
            let mut best = None::<(isize, u8)>;
            // Equal furthest-reaching states use X, I, D in that stable order.
            for (previous_k, op) in [(k, b'X'), (k - 1, b'I'), (k + 1, b'D')] {
                if let Some(q) = previous.at(score - 1, previous_k) {
                    let t = q - previous_k;
                    let legal = match op {
                        b'X' => q < m && t < n,
                        b'I' => q < m,
                        _ => t < n,
                    };
                    if !legal {
                        continue;
                    }
                    let next = q + isize::from(op != b'D');
                    if best.is_none_or(|(position, _)| next > position) {
                        best = Some((next, op));
                    }
                }
            }
            if let Some((q, op)) = best {
                let index = (k + score as isize) as usize;
                layer.query[index] = extend(q, k)?;
                layer.transition[index] = op;
            }
        }
        layers.push(layer);
    }
    let mut q = m;
    let mut k = end_diagonal;
    let mut reversed = Vec::new();
    for s in (1..=score).rev() {
        let op = layers[s].transition[(k + s as isize) as usize];
        let previous_k = match op {
            b'X' => k,
            b'I' => k - 1,
            b'D' => k + 1,
            _ => return Err("global_wavefront_traceback_invariant".to_string()),
        };
        let previous_q = layers[s - 1]
            .at(s - 1, previous_k)
            .ok_or("global_wavefront_traceback_invariant")?;
        let after_edit = previous_q + isize::from(op != b'D');
        push_op(
            &mut reversed,
            usize::try_from(q - after_edit).map_err(|_| "global_wavefront_negative_match_span")?,
            b'=',
        );
        push_op(&mut reversed, 1, op);
        q = previous_q;
        k = previous_k;
    }
    push_op(&mut reversed, q as usize, b'=');
    reversed.reverse();
    Ok(GlobalAlignment {
        cigar: normalized_cigar(source, target, &reversed)?,
        edit_distance: score,
        explored_cells,
    })
}

fn global_candidates_from_paf(
    source: &str,
    target: &str,
    circular: bool,
    paf: &str,
) -> Result<Vec<AlignmentCandidate>, String> {
    let reference_len = target
        .len()
        .checked_mul(if circular { 2 } else { 1 })
        .ok_or("global_reference_length_overflow")?;
    let mut offsets = BTreeSet::new();
    let mut reasons = BTreeSet::new();
    for line in paf.lines().filter(|line| !line.trim().is_empty()) {
        let anchor = match parse_anchor(line, source.len(), reference_len) {
            Ok(anchor) => anchor,
            Err(reason) => {
                reasons.insert(reason);
                continue;
            }
        };
        let core_target = periodic_slice(
            target.as_bytes(),
            anchor.target_start,
            anchor.target_end - anchor.target_start,
        );
        // The fallback does not excuse false '='/'X' evidence in the original anchor.
        normalized_cigar(
            &source.as_bytes()[anchor.query_start..anchor.query_end],
            &core_target,
            &anchor.cigar,
        )?;
        if circular {
            for (t, q) in [
                (anchor.target_start, anchor.query_start),
                (anchor.target_end, anchor.query_end),
            ] {
                offsets.insert((t % target.len() + target.len() - q % target.len()) % target.len());
            }
        } else {
            offsets.insert(0);
        }
    }
    if offsets.is_empty() {
        return Err(
            "no_complete_candidate: no_valid_forward_anchor_for_global_alignment".to_string(),
        );
    }
    if offsets.len() > MAX_CANDIDATES {
        return Err("no_complete_candidate: implementation_limit:global_frame_limit".to_string());
    }
    let mut candidates = Vec::new();
    for offset in offsets {
        let frame = rotate(target, offset);
        match complete_global_cigar(source.as_bytes(), frame.as_bytes()) {
            Ok(alignment) => candidates.push(AlignmentCandidate {
                target_in_path_frame: frame, frame_offset: offset, cigar: alignment.cigar,
                terminal_target_overlap: 0,
                method: format!("complete_global_wavefront;canonical_optimal_traceback;edit_distance={};explored_cells={}", alignment.edit_distance, alignment.explored_cells),
            }),
            Err(reason) => { reasons.insert(format!("frame_{offset}:{reason}")); }
        }
    }
    if candidates.is_empty() {
        return Err(format!(
            "no_complete_candidate: global_alignment_implementation_limit:{}",
            reasons.into_iter().collect::<Vec<_>>().join(";")
        ));
    }
    if !reasons.is_empty() {
        let limitations = reasons.into_iter().collect::<Vec<_>>().join(",");
        for candidate in &mut candidates {
            candidate
                .method
                .push_str(&format!(";other_frames_not_completed={limitations}"));
        }
    }
    Ok(candidates)
}

fn parse_anchor(line: &str, source_len: usize, target_len: usize) -> Result<Anchor, String> {
    let f: Vec<_> = line.split('\t').collect();
    if f.len() < 12 || f[0] != "path_source" || f[5] != "path_target" || f[4] != "+" {
        return Err("not_a_forward_fixed_path_alignment".to_string());
    }
    let number = |index: usize| -> Result<usize, String> {
        f[index]
            .parse::<usize>()
            .map_err(|_| "invalid_path_paf_coordinate".to_string())
    };
    if number(1)? != source_len || number(6)? != target_len {
        return Err("path_paf_sequence_lengths_disagree".to_string());
    }
    let (qs, qe, ts, te) = (number(2)?, number(3)?, number(7)?, number(8)?);
    if qs >= qe || qe > source_len || ts >= te || te > target_len {
        return Err("path_paf_coordinates_out_of_range".to_string());
    }
    let cigar = f
        .iter()
        .find_map(|v| v.strip_prefix("cg:Z:"))
        .ok_or("path_paf_lacks_cigar")?;
    let cigar = parse_cigar(cigar)?;
    if consumed(&cigar)? != (qe - qs, te - ts) {
        return Err("path_paf_cigar_span_disagrees".to_string());
    }
    Ok(Anchor {
        query_start: qs,
        query_end: qe,
        target_start: ts,
        target_end: te,
        cigar,
        mapq: number(11)?,
        primary: f.iter().any(|field| *field == "tp:A:P"),
    })
}

struct Matrix<'a> {
    query: &'a [u8],
    target: &'a [u8],
    cols: usize,
    scores: Vec<u32>,
}
impl<'a> Matrix<'a> {
    fn new(query: &'a [u8], target: &'a [u8]) -> Result<Self, String> {
        if query.len() > MAX_TERMINAL_QUERY_BP || target.len() > MAX_TERMINAL_TARGET_BP {
            return Err("terminal_extension_length_limit".to_string());
        }
        let cells = (query.len() + 1)
            .checked_mul(target.len() + 1)
            .ok_or("terminal_extension_size_overflow")?;
        if cells > MAX_DP_CELLS {
            return Err("terminal_extension_memory_limit".to_string());
        }
        let cols = target.len() + 1;
        let mut scores = vec![0u32; cells];
        for i in 0..=query.len() {
            scores[i * cols] = i as u32;
        }
        for j in 0..=target.len() {
            scores[j] = j as u32;
        }
        for i in 1..=query.len() {
            for j in 1..=target.len() {
                scores[i * cols + j] = (scores[(i - 1) * cols + j - 1]
                    + u32::from(query[i - 1] != target[j - 1]))
                .min(scores[(i - 1) * cols + j] + 1)
                .min(scores[i * cols + j - 1] + 1);
            }
        }
        Ok(Self {
            query,
            target,
            cols,
            scores,
        })
    }
    fn score(&self, i: usize, j: usize) -> u32 {
        self.scores[i * self.cols + j]
    }
    fn tracebacks(&self, target_length: usize) -> Result<Vec<Vec<Op>>, String> {
        let mut pending = vec![(self.query.len(), target_length, Vec::<u8>::new())];
        let mut results = Vec::new();
        let mut visited = 0;
        while let Some((i, j, reverse_ops)) = pending.pop() {
            visited += 1;
            if visited > MAX_TRACEBACK_STATES {
                return Err("terminal_extension_traceback_limit".to_string());
            }
            if i == 0 && j == 0 {
                let mut ops = Vec::new();
                for &code in reverse_ops.iter().rev() {
                    push_op(&mut ops, 1, code);
                }
                results.push(ops);
                if results.len() > MAX_CANDIDATES {
                    return Err("terminal_extension_has_too_many_equal_alignments".to_string());
                }
                continue;
            }
            let best = self.score(i, j);
            let mut choices = Vec::new();
            if i > 0 && j > 0 {
                let mismatch = self.query[i - 1] != self.target[j - 1];
                if self.score(i - 1, j - 1) + u32::from(mismatch) == best {
                    choices.push((i - 1, j - 1, if mismatch { b'X' } else { b'=' }));
                }
            }
            if i > 0 && self.score(i - 1, j) + 1 == best {
                choices.push((i - 1, j, b'I'));
            }
            if j > 0 && self.score(i, j - 1) + 1 == best {
                choices.push((i, j - 1, b'D'));
            }
            for (ni, nj, code) in choices {
                let mut next_ops = reverse_ops.clone();
                next_ops.push(code);
                pending.push((ni, nj, next_ops));
            }
            if pending.len() > MAX_CANDIDATES {
                return Err("terminal_extension_pending_candidate_limit".to_string());
            }
        }
        Ok(results)
    }
}

fn unique_circular_flank(target: &str, flank: &[u8]) -> bool {
    let doubled = format!("{target}{target}");
    doubled
        .as_bytes()
        .windows(flank.len())
        .take(target.len())
        .filter(|window| *window == flank)
        .take(2)
        .count()
        == 1
}

fn circular_seam_gap_candidate(
    source: &str,
    target: &str,
    anchor: &Anchor,
    slack: usize,
) -> Option<AlignmentCandidate> {
    let target_len = target.len();
    if slack == 0
        || slack > MAX_TERMINAL_TARGET_BP
        || !anchor.primary
        || anchor.mapq < 20
        || anchor.query_start != 0
        || anchor.query_end != source.len()
        || anchor.target_start >= target_len
        || anchor.target_end <= target_len
        || anchor.target_end > target_len.checked_mul(2)?
        || !anchor.cigar.first().is_some_and(|op| op.code == b'=' && op.count >= UNIQUE_PROJECTION_FLANK_BP)
        || !anchor.cigar.last().is_some_and(|op| op.code == b'=' && op.count >= UNIQUE_PROJECTION_FLANK_BP)
    {
        return None;
    }
    let offset = anchor.target_start;
    let target_frame = rotate(target, offset);
    let cut = target_len - slack;
    if cut < UNIQUE_PROJECTION_FLANK_BP || source.len() < UNIQUE_PROJECTION_FLANK_BP * 2 {
        return None;
    }
    let left = &target_frame.as_bytes()[cut - UNIQUE_PROJECTION_FLANK_BP..cut];
    let right = &target_frame.as_bytes()[..UNIQUE_PROJECTION_FLANK_BP];
    if !unique_circular_flank(target, left) || !unique_circular_flank(target, right) {
        return None;
    }
    let mut ops = anchor.cigar.clone();
    push_op(&mut ops, slack, b'D');
    let cigar = normalized_cigar(source.as_bytes(), target_frame.as_bytes(), &ops).ok()?;
    Some(AlignmentCandidate {
        target_in_path_frame: target_frame,
        frame_offset: offset,
        cigar,
        method: "minimap2_primary_complete_path_circular_seam_gap".to_string(),
        terminal_target_overlap: 0,
    })
}

/// A full-query anchor may spell one circle plus its repeated closing prefix.
/// Keep that complete witness; the graph projection computes its overlap later.
fn circular_seam_overlap_candidate(
    source: &str,
    target: &str,
    anchor: &Anchor,
) -> Option<AlignmentCandidate> {
    let span = anchor.target_end.checked_sub(anchor.target_start)?;
    let overlap = span.checked_sub(target.len())?;
    if overlap == 0
        || overlap > MAX_TERMINAL_TARGET_BP
        || overlap >= target.len()
        || overlap.checked_mul(2)? >= source.len()
        || !anchor.primary
        || anchor.mapq < 20
        || anchor.query_start != 0
        || anchor.query_end != source.len()
        || anchor.target_start >= target.len()
        || anchor.target_end > target.len().checked_mul(2)?
        || !anchor.cigar.first().is_some_and(|op| {
            op.code == b'=' && op.count >= overlap.max(UNIQUE_PROJECTION_FLANK_BP)
        })
        || !anchor.cigar.last().is_some_and(|op| {
            op.code == b'=' && op.count >= overlap.max(UNIQUE_PROJECTION_FLANK_BP)
        })
        || source.as_bytes()[..overlap] != source.as_bytes()[source.len() - overlap..]
    {
        return None;
    }
    let frame = rotate(target, anchor.target_start);
    // A short seam can recur internally. Uniqueness concerns the entire circle,
    // not an arbitrary flank that may itself belong to a genuine repeat node.
    if !unique_circular_flank(target, frame.as_bytes())
        || frame.as_bytes()[..overlap] != source.as_bytes()[..overlap]
    {
        return None;
    }
    let unfolded = format!("{frame}{}", &frame[..overlap]);
    let cigar = normalized_cigar(source.as_bytes(), unfolded.as_bytes(), &anchor.cigar).ok()?;
    Some(AlignmentCandidate {
        target_in_path_frame: frame,
        frame_offset: anchor.target_start,
        cigar,
        method: format!(
            "minimap2_primary_complete_path_circular_seam_overlap;terminal_target_overlap={overlap};unfolded_target_length={span};query_clipped=false"
        ),
        terminal_target_overlap: overlap,
    })
}

fn candidates_from_anchor(
    source: &str,
    target: &str,
    circular: bool,
    anchor: &Anchor,
) -> Result<Vec<AlignmentCandidate>, String> {
    let left_query = &source.as_bytes()[..anchor.query_start];
    let right_query = &source.as_bytes()[anchor.query_end..];
    let core_target_length = anchor.target_end - anchor.target_start;
    if left_query.len() > MAX_TERMINAL_QUERY_BP || right_query.len() > MAX_TERMINAL_QUERY_BP {
        return Err("path_anchor_leaves_long_source_ends".to_string());
    }
    if core_target_length > target.len() {
        if circular {
            if let Some(candidate) = circular_seam_overlap_candidate(source, target, anchor) {
                return Ok(vec![candidate]);
            }
        }
        return Err("path_anchor_consumes_more_than_one_target_copy".to_string());
    }
    let slack = target.len() - core_target_length;
    if slack > MAX_TERMINAL_TARGET_BP {
        return Err("path_anchor_leaves_long_target_ends".to_string());
    }
    if circular {
        if let Some(candidate) = circular_seam_gap_candidate(source, target, anchor, slack) {
            return Ok(vec![candidate]);
        }
    }
    let left_query_reverse: Vec<_> = left_query.iter().rev().copied().collect();
    let (left_target_reverse, right_target, allowed_splits): (Vec<u8>, Vec<u8>, Vec<usize>) =
        if circular {
            let left_start = (anchor.target_start % target.len() + target.len()
                - slack % target.len())
                % target.len();
            let left = periodic_slice(target.as_bytes(), left_start, slack);
            (
                left.into_iter().rev().collect(),
                periodic_slice(target.as_bytes(), anchor.target_end % target.len(), slack),
                (0..=slack).collect(),
            )
        } else {
            (
                target.as_bytes()[..anchor.target_start]
                    .iter()
                    .rev()
                    .copied()
                    .collect(),
                target.as_bytes()[anchor.target_end..].to_vec(),
                vec![anchor.target_start],
            )
        };
    let prefix = Matrix::new(&left_query_reverse, &left_target_reverse)?;
    let suffix = Matrix::new(right_query, &right_target)?;
    let best = allowed_splits
        .iter()
        .map(|&left| {
            let right = slack - left;
            prefix.score(left_query.len(), left) + suffix.score(right_query.len(), right)
        })
        .min()
        .ok_or("no_terminal_extension_split")?;
    let mut candidates = Vec::new();
    for left in allowed_splits {
        let right = slack - left;
        if prefix.score(left_query.len(), left) + suffix.score(right_query.len(), right) != best {
            continue;
        }
        let offset = if circular {
            (anchor.target_start % target.len() + target.len() - left % target.len()) % target.len()
        } else {
            0
        };
        let target_frame = rotate(target, offset);
        let prefixes = prefix.tracebacks(left)?;
        let suffixes = suffix.tracebacks(right)?;
        for pre_reverse in &prefixes {
            let mut pre = pre_reverse.clone();
            pre.reverse();
            for post in &suffixes {
                let mut full = Vec::new();
                append_ops(&mut full, &pre);
                append_ops(&mut full, &anchor.cigar);
                append_ops(&mut full, post);
                let cigar = normalized_cigar(source.as_bytes(), target_frame.as_bytes(), &full)?;
                candidates.push(AlignmentCandidate {
                    target_in_path_frame: target_frame.clone(),
                    frame_offset: offset,
                    cigar,
                    method: if left_query.is_empty() && right_query.is_empty() && slack == 0 {
                        "minimap2_complete_path"
                    } else {
                        "minimap2_with_bounded_global_end_extension"
                    }
                    .to_string(),
                    terminal_target_overlap: 0,
                });
                if candidates.len() > MAX_CANDIDATES {
                    return Err("path_alignment_candidate_limit".to_string());
                }
            }
        }
    }
    Ok(candidates)
}

fn candidates_from_paf(
    source: &str,
    target: &str,
    circular: bool,
    paf: &str,
) -> Result<Vec<AlignmentCandidate>, String> {
    let reference_len = if circular {
        target
            .len()
            .checked_mul(2)
            .ok_or("doubled_target_length_overflow")?
    } else {
        target.len()
    };
    let mut candidates = Vec::new();
    let mut seen = BTreeSet::new();
    let mut reasons = BTreeSet::new();
    for line in paf.lines().filter(|line| !line.trim().is_empty()) {
        let anchor = match parse_anchor(line, source.len(), reference_len) {
            Ok(anchor) => anchor,
            Err(reason) => {
                reasons.insert(reason);
                continue;
            }
        };
        match candidates_from_anchor(source, target, circular, &anchor) {
            Ok(found) => {
                for candidate in found {
                    if seen.insert((
                        candidate.target_in_path_frame.clone(),
                        candidate.cigar.clone(),
                    )) {
                        candidates.push(candidate);
                    }
                    if candidates.len() > MAX_CANDIDATES {
                        return Err(
                            "no_complete_candidate: resource_exhausted:path_alignment_candidate_limit".to_string()
                        );
                    }
                }
            }
            Err(reason) => {
                // A bounded ambiguity search that could not finish must not be
                // hidden by accepting a different anchor's convenient candidate.
                if reason.contains("_limit") || reason.contains("too_many_equal_alignments") {
                    return Err(format!(
                        "no_complete_candidate: resource_exhausted:candidate_enumeration_incomplete:{reason}"
                    ));
                }
                reasons.insert(reason);
            }
        }
    }
    if candidates.iter().any(|candidate| candidate.terminal_target_overlap > 0)
        && candidates.len() != 1
    {
        return Err("circular_seam_has_multiple_complete_alignment_witnesses".to_string());
    }
    if candidates.is_empty() {
        let detail = if reasons.is_empty() {
            "minimap2_returned_no_forward_path_hit".to_string()
        } else {
            reasons.into_iter().collect::<Vec<_>>().join(";")
        };
        return Err(format!("no_complete_candidate: {detail}"));
    }
    if !reasons.is_empty() {
        let limitations = reasons.into_iter().collect::<Vec<_>>().join(",");
        for candidate in &mut candidates {
            candidate.method.push_str("; other_anchors_not_completed=");
            candidate.method.push_str(&limitations);
        }
    }
    Ok(candidates)
}

#[cfg(test)]
mod tests_path_alignment {
    use super::*;
    fn paf(
        q_len: usize,
        qs: usize,
        qe: usize,
        t_len: usize,
        ts: usize,
        te: usize,
        cigar: &str,
    ) -> String {
        format!("path_source\t{q_len}\t{qs}\t{qe}\t+\tpath_target\t{t_len}\t{ts}\t{te}\t10\t10\t60\tcg:Z:{cigar}")
    }

    fn dense_edit_distance(query: &[u8], target: &[u8]) -> usize {
        let mut previous = (0..=target.len()).collect::<Vec<_>>();
        for (i, q) in query.iter().enumerate() {
            let mut next = vec![i + 1; target.len() + 1];
            for (j, t) in target.iter().enumerate() {
                next[j + 1] = (previous[j] + usize::from(q != t))
                    .min(previous[j + 1] + 1)
                    .min(next[j] + 1);
            }
            previous = next;
        }
        previous[target.len()]
    }

    #[test]
    fn global_wavefront_matches_independent_dense_distance_and_full_consumption() {
        let mut sequences = vec![String::new()];
        for n in 1..=5 {
            for bits in 0..(1 << n) {
                sequences.push(
                    (0..n)
                        .map(|i| if (bits >> i) & 1 == 0 { 'A' } else { 'C' })
                        .collect(),
                );
            }
        }
        for query in &sequences {
            for target in &sequences {
                let alignment = complete_global_cigar(query.as_bytes(), target.as_bytes()).unwrap();
                assert_eq!(
                    alignment.edit_distance,
                    dense_edit_distance(query.as_bytes(), target.as_bytes())
                );
                if !query.is_empty() || !target.is_empty() {
                    let ops = parse_cigar(&alignment.cigar).unwrap();
                    assert_eq!(consumed(&ops).unwrap(), (query.len(), target.len()));
                    assert_eq!(
                        normalized_cigar(query.as_bytes(), target.as_bytes(), &ops).unwrap(),
                        alignment.cigar
                    );
                }
            }
        }
    }

    #[test]
    fn a_complete_global_alignment_can_recover_a_long_clipped_source_end() {
        let source = format!("{}CGTA", "A".repeat(MAX_TERMINAL_QUERY_BP + 20));
        let mut target = source.clone();
        target.push('C');
        let hit = paf(source.len(), 0, 10, target.len(), 0, 10, "10=");
        let old_error = candidates_from_paf(&source, &target, false, &hit).unwrap_err();
        assert!(can_try_complete_global_alignment(&old_error));
        let result = global_candidates_from_paf(&source, &target, false, &hit).unwrap();
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].cigar, format!("{}=1D", source.len()));
    }

    #[test]
    fn a_one_base_anchor_crossing_the_circle_is_realigned_to_exactly_one_target() {
        let hit = paf(5, 0, 5, 8, 0, 5, "5=");
        assert!(candidates_from_paf("ACGTA", "ACGT", true, &hit).is_err());
        let result = global_candidates_from_paf("ACGTA", "ACGT", true, &hit).unwrap();
        assert_eq!(result[0].cigar, "4=1I");
        assert_eq!(
            consumed(&parse_cigar(&result[0].cigar).unwrap()).unwrap(),
            (5, 4)
        );
    }

    #[test]
    fn canonical_global_traceback_does_not_enumerate_exponentially_many_equal_gaps() {
        let source = format!("{}C", "A".repeat(12));
        let target = format!("{}C", "A".repeat(24));
        let hit = paf(13, 12, 13, 25, 24, 25, "1=");
        let error = candidates_from_paf(&source, &target, false, &hit).unwrap_err();
        assert!(error.contains("candidate_enumeration_incomplete"));
        let result = global_candidates_from_paf(&source, &target, false, &hit).unwrap();
        assert_eq!(result.len(), 1);
        assert!(result[0]
            .method
            .contains("canonical_optimal_traceback;edit_distance=12"));
        assert_eq!(
            consumed(&parse_cigar(&result[0].cigar).unwrap()).unwrap(),
            (13, 25)
        );
    }

    #[test]
    fn global_fallback_keeps_distinct_anchor_frames_and_rejects_bad_equality_evidence() {
        let hits = format!(
            "{}\n{}",
            paf(5, 0, 2, 12, 0, 2, "2="),
            paf(5, 2, 5, 12, 3, 6, "3=")
        );
        let result = global_candidates_from_paf("AACGT", "AAACGT", true, &hits).unwrap();
        assert_eq!(
            result
                .iter()
                .map(|c| c.frame_offset)
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([0, 1])
        );
        assert!(result
            .iter()
            .all(|c| consumed(&parse_cigar(&c.cigar).unwrap()).unwrap() == (5, 6)));
        let bad = paf(4, 0, 4, 4, 0, 4, "4X");
        assert!(global_candidates_from_paf("AAAA", "AAAA", false, &bad)
            .unwrap_err()
            .contains("alignment_mismatch_operator_conflicts_with_sequence"));
        assert!(!can_try_complete_global_alignment(
            "no_complete_candidate: alignment_equality_conflicts_with_sequence"
        ));
    }

    #[test]
    fn global_alignment_resource_limits_remain_explicit() {
        assert_eq!(
            complete_global_cigar_with_limits(b"AAA", b"CCC", 3, 100)
                .err()
                .unwrap(),
            "global_wavefront_memory_limit"
        );
        assert_eq!(
            complete_global_cigar_with_limits(b"AAAA", b"AAAA", 100, 1)
                .err()
                .unwrap(),
            "global_wavefront_work_limit"
        );
    }
    #[test]
    fn exact_shortcuts_do_not_call_a_tool_or_write_files() {
        let p = Path::new("/not/a/minimap2");
        let work = Path::new("/not/a/workdir");
        let c = align_path("ACGTAA", "GTAAAC", true, false, p, 1, work).unwrap();
        assert_eq!(c[0].target_in_path_frame, "ACGTAA");
        assert_eq!(c[0].frame_offset, 4);
        let c = align_path("AACC", "GGTT", false, true, p, 1, work).unwrap();
        assert_eq!(c[0].cigar, "4=");
    }
    #[test]
    fn complete_paf_uses_source_query_indel_semantics() {
        let c = candidates_from_paf(
            "AACGT",
            "AATCGTT",
            false,
            &paf(5, 0, 5, 7, 0, 7, "2=1D3=1D"),
        )
        .unwrap();
        assert_eq!(c[0].cigar, "2=1D3=1D");
    }
    #[test]
    fn fills_both_linear_ends_without_dropping_bases() {
        let c =
            candidates_from_paf("AACGTT", "AATCGTTA", false, &paf(6, 2, 5, 8, 3, 6, "3=")).unwrap();
        assert!(c
            .iter()
            .all(|x| consumed(&parse_cigar(&x.cigar).unwrap()).unwrap() == (6, 8)));
        assert!(c.iter().any(|x| x.cigar == "2=1D4=1D"));
    }
    #[test]
    fn restores_circular_prefix_before_anchor_at_zero() {
        let c = candidates_from_paf("TAACG", "AACGT", true, &paf(5, 1, 5, 10, 0, 4, "4=")).unwrap();
        assert!(c.iter().any(|x| x.frame_offset == 4 && x.cigar == "5="));
    }
    #[test]
    fn retains_tied_terminal_alignments() {
        let c = candidates_from_paf("AACG", "AAACG", false, &paf(4, 2, 4, 5, 3, 5, "2=")).unwrap();
        let set: BTreeSet<_> = c.iter().map(|x| x.cigar.as_str()).collect();
        assert_eq!(set, BTreeSet::from(["1D4=", "1=1D3=", "2=1D2="]));
    }
    #[test]
    fn rejects_partial_hit_with_a_large_unaligned_tail() {
        let source = "A".repeat(MAX_TERMINAL_QUERY_BP + 20);
        let error = candidates_from_paf(
            &source,
            &source,
            false,
            &paf(source.len(), 0, 10, source.len(), 0, 10, "10="),
        )
        .unwrap_err();
        assert!(error.starts_with("no_complete_candidate:"));
        assert!(error.contains("long_source_ends"));
    }
    #[test]
    fn cannot_hide_an_unfinished_equal_candidate_search_with_another_anchor() {
        let source = format!("{}C", "A".repeat(12));
        let target = format!("{}C", "A".repeat(24));
        let complete = paf(13, 0, 13, 25, 0, 25, "12D13=");
        let clipped = paf(13, 12, 13, 25, 24, 25, "1=");
        let error = candidates_from_paf(&source, &target, false, &format!("{complete}\n{clipped}"))
            .unwrap_err();
        assert!(error.contains("candidate_enumeration_incomplete"));
    }

    #[test]
    fn rejects_explicit_mismatch_when_the_bases_are_equal() {
        let error =
            candidates_from_paf("AAAA", "AAAA", false, &paf(4, 0, 4, 4, 0, 4, "4X")).unwrap_err();
        assert!(error.contains("alignment_mismatch_operator_conflicts_with_sequence"));
    }

    #[test]
    fn reports_skipped_long_partial_hits_beside_a_complete_candidate() {
        let source = "A".repeat(MAX_TERMINAL_QUERY_BP + 20);
        let n = source.len();
        let complete = paf(n, 0, n, n, 0, n, &format!("{n}="));
        let partial = paf(n, 0, 10, n, 0, 10, "10=");
        let result =
            candidates_from_paf(&source, &source, false, &format!("{complete}\n{partial}"))
                .unwrap();
        assert_eq!(result.len(), 1);
        assert!(result[0]
            .method
            .contains("other_anchors_not_completed=path_anchor_leaves_long_source_ends"));
    }

    #[test]
    fn rejects_reverse_hits_and_more_than_one_target_copy() {
        let reverse = paf(4, 0, 4, 8, 0, 4, "4=").replace("\t+\t", "\t-\t");
        assert!(candidates_from_paf("ACGT", "ACGT", true, &reverse).is_err());
        assert!(
            candidates_from_paf("ACGTACGT", "ACGT", true, &paf(8, 0, 8, 8, 0, 8, "8=")).is_err()
        );
    }
}



// =====================================================================
// Inlined module: rebuild_path_sequence.rs (walk projection through alignment)
// =====================================================================
// Project a fixed graph walk through a complete alignment to an authoritative sequence.
//
// This module never selects a walk, changes node identity, or splits a shared node.
// The caller supplies the target in the original walk's direction and coordinate frame.
// Alignment CIGAR uses the source spell as query and the target as reference: I consumes
// source only, D consumes target only. A successful plan spells the target exactly and
// checks the sequence overlap of every supplied graph edge, including unused edges.


#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PathStep {
    pub node: String,
    pub reverse: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PathEdge {
    pub from: String,
    pub from_reverse: bool,
    pub to: String,
    pub to_reverse: bool,
    /// Number of matching bases in a simple nM/n= overlap. The caller must reject
    /// complex or unknown CIGARs rather than converting them to zero.
    pub overlap: usize,
}

#[derive(Debug, Clone)]
struct SourceOccurrence {
    step: PathStep,
    start: isize,
    end: isize,
    incoming_edges: Vec<usize>,
}

/// A source-graph edge whose declared overlap end sequences do not match exactly.
/// The recorded walk still uses the declared overlap to spell the source query; the
/// authoritative target projection must afterwards make EVERY output edge exact. A
/// recorded conflict is never repaired by deleting source bases on its own.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceOverlapConflict {
    pub edge_index: usize,
    pub from: String,
    pub to: String,
    pub declared_overlap: usize,
    pub mismatch_count: usize,
    /// Zero-based positions inside the declared overlap window.
    pub mismatch_positions: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct PreparedPath {
    /// For a circle, each occurrence contributes after its incoming overlap, including
    /// the first occurrence. Thus the closing overlap is represented exactly once.
    pub source_spell: String,
    pub circular: bool,
    edges: Vec<PathEdge>,
    occurrences: Vec<SourceOccurrence>,
    /// Source overlap conflicts that were tolerated to build this source spell. Only
    /// the corrected (lenient) preparation populates this; strict callers reject first.
    pub source_overlap_conflicts: Vec<SourceOverlapConflict>,
}

impl PreparedPath {
    /// Occurrence intervals use the same zero-based half-open convention as the output
    /// copies, but refer to source_spell. Circular intervals are deliberately unfolded.
    pub fn source_copies(&self) -> Vec<CopyPlacement> {
        self.occurrences
            .iter()
            .map(|copy| CopyPlacement {
                node: copy.step.node.clone(),
                reverse: copy.step.reverse,
                start: copy.start,
                end: copy.end,
            })
            .collect()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CopyPlacement {
    pub node: String,
    pub reverse: bool,
    /// Zero-based half-open positions in the supplied target frame. Circular copies
    /// may start below zero; periodic extraction, not endpoint-only modulo, is required.
    pub start: isize,
    pub end: isize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdgeOverlap {
    pub edge_index: usize,
    pub old_overlap: usize,
    pub new_overlap: usize,
    pub supported_by_walk_alignment: bool,
}

#[derive(Debug, Clone)]
pub struct SharedCoreTrim {
    pub node: String,
    pub source_length: usize,
    /// Zero-based, half-open coordinates in the original node's forward sequence.
    pub source_start: usize,
    pub source_end: usize,
    pub new_length: usize,
}

#[derive(Debug, Clone)]
pub struct PathSequencePlan {
    pub sequences: BTreeMap<String, String>,
    pub copies: Vec<CopyPlacement>,
    pub edge_overlaps: Vec<EdgeOverlap>,
    pub explored_assignments: usize,
    /// Complete alignment actually used for the accepted boundary assignment.
    pub alignment_cigar: String,
    /// Includes the supplied alignment and any equivalent indel placements tried.
    pub alignment_candidates_examined: usize,
    /// Rotation of the accepted graph frame relative to the alignment target frame.
    /// Nonzero values are permitted only for an explicitly circular path.
    pub target_frame_shift: isize,
    pub shared_core_trims: Vec<SharedCoreTrim>,
}

fn reverse_complement_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'M' => b'K',
            b'K' => b'M',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            _ => *base,
        })
        .collect()
}

fn dna(sequence: &str) -> Result<Vec<u8>, String> {
    let sequence = sequence.to_ascii_uppercase().into_bytes();
    if sequence.is_empty() || sequence.iter().any(|b| !b"ACGTRYSWKMBDHVN".contains(b)) {
        return Err("sequence_is_empty_or_contains_non_iupac_dna".to_string());
    }
    Ok(sequence)
}

fn oriented_bytes(sequence: &[u8], reverse: bool) -> Vec<u8> {
    if reverse {
        reverse_complement_bytes(sequence)
    } else {
        sequence.to_vec()
    }
}

fn matches_edge(edge: &PathEdge, left: &PathStep, right: &PathStep) -> bool {
    (edge.from == left.node
        && edge.from_reverse == left.reverse
        && edge.to == right.node
        && edge.to_reverse == right.reverse)
        || (edge.to == left.node
            && edge.to_reverse != left.reverse
            && edge.from == right.node
            && edge.from_reverse != right.reverse)
}

/// Prepare only an already chosen, legal, complete walk. Node copies in the walk use
/// the original shared node ID; the caller resolves historical copy-instance aliases.
/// Strict source preparation: every declared source overlap must be an exact terminal
/// match. Use this for the node-projection generator and all structural checks.
pub fn prepare_path(
    nodes: &BTreeMap<String, String>,
    edges: &[PathEdge],
    walk: &[PathStep],
    circular: bool,
) -> Result<PreparedPath, String> {
    prepare_path_inner(nodes, edges, walk, circular, false)
}

/// Lenient preparation used only by the whole-path alignment generator when the source
/// graph carries documented overlap end-sequence errors that polish corrected. Conflicts
/// are recorded on the returned path; the target projection in `evaluate` still requires
/// every output edge (walked and unused) to be exactly compatible and the full replay to
/// equal the authoritative FASTA, so a recorded conflict can never pass on its own.
pub fn prepare_path_corrected(
    nodes: &BTreeMap<String, String>,
    edges: &[PathEdge],
    walk: &[PathStep],
    circular: bool,
) -> Result<PreparedPath, String> {
    prepare_path_inner(nodes, edges, walk, circular, true)
}

fn prepare_path_inner(
    nodes: &BTreeMap<String, String>,
    edges: &[PathEdge],
    walk: &[PathStep],
    circular: bool,
    tolerate_overlap_conflict: bool,
) -> Result<PreparedPath, String> {
    if walk.is_empty() {
        return Err("fixed_walk_is_empty".to_string());
    }
    let mut overlap_conflicts: Vec<SourceOverlapConflict> = Vec::new();
    let nodes = nodes
        .iter()
        .map(|(name, seq)| Ok((name.clone(), dna(seq)?)))
        .collect::<Result<BTreeMap<_, _>, String>>()?;
    let used = walk
        .iter()
        .map(|step| step.node.clone())
        .collect::<BTreeSet<_>>();
    if used != nodes.keys().cloned().collect() {
        return Err("fixed_walk_does_not_cover_exactly_the_component_nodes".to_string());
    }
    for edge in edges {
        if !nodes.contains_key(&edge.from) || !nodes.contains_key(&edge.to) {
            return Err("edge_references_node_outside_component".to_string());
        }
    }
    let mut source = Vec::new();
    let mut occurrences = Vec::new();
    for (index, step) in walk.iter().enumerate() {
        let incoming_edges = if index > 0 || circular {
            let previous = &walk[(index + walk.len() - 1) % walk.len()];
            let found = edges
                .iter()
                .enumerate()
                .filter_map(|(i, edge)| matches_edge(edge, previous, step).then_some(i))
                .collect::<Vec<_>>();
            if found.is_empty() {
                return Err(format!(
                    "fixed_walk_has_no_oriented_edge:{}:{}",
                    previous.node, step.node
                ));
            }
            let overlaps = found
                .iter()
                .map(|i| edges[*i].overlap)
                .collect::<BTreeSet<_>>();
            if overlaps.len() != 1 {
                return Err(format!(
                    "fixed_walk_has_ambiguous_edge_overlap:{}:{}",
                    previous.node, step.node
                ));
            }
            let n = edges[found[0]].overlap;
            if n > nodes[&previous.node].len() || n > nodes[&step.node].len() {
                return Err("source_edge_overlap_exceeds_node_length".to_string());
            }
            let left = oriented_bytes(&nodes[&previous.node], previous.reverse);
            let right = oriented_bytes(&nodes[&step.node], step.reverse);
            if left[left.len() - n..] != right[..n] {
                let window_left = &left[left.len() - n..];
                let window_right = &right[..n];
                let positions = window_left
                    .iter()
                    .zip(window_right.iter())
                    .enumerate()
                    .filter_map(|(pos, (a, b))| (a != b).then_some(pos))
                    .collect::<Vec<_>>();
                if !tolerate_overlap_conflict || (index == 0 && circular) {
                    if index == 0 && circular {
                        return Err(
                            "source_closing_overlap_is_not_an_exact_terminal_match".to_string()
                        );
                    }
                    return Err(format!(
                        "source_edge_overlap_sequence_conflict:{}:{}:{}",
                        previous.node, step.node, n
                    ));
                }
                for edge_index in &found {
                    overlap_conflicts.push(SourceOverlapConflict {
                        edge_index: *edge_index,
                        from: previous.node.clone(),
                        to: step.node.clone(),
                        declared_overlap: n,
                        mismatch_count: positions.len(),
                        mismatch_positions: positions.clone(),
                    });
                }
            }
            found
        } else {
            Vec::new()
        };
        let overlap = incoming_edges.first().map_or(0, |i| edges[*i].overlap);
        let sequence = oriented_bytes(&nodes[&step.node], step.reverse);
        let start = source.len() as isize - overlap as isize;
        source.extend_from_slice(&sequence[overlap..]);
        occurrences.push(SourceOccurrence {
            step: step.clone(),
            start,
            end: source.len() as isize,
            incoming_edges,
        });
    }
    if source.is_empty() {
        return Err("source_walk_has_no_sequence_after_overlaps".to_string());
    }
    if circular {
        let closing_overlap = edges[occurrences[0].incoming_edges[0]].overlap;
        if closing_overlap > 0 {
            let first = oriented_bytes(&nodes[&walk[0].node], walk[0].reverse);
            let mut untrimmed = first[..closing_overlap].to_vec();
            untrimmed.extend_from_slice(&source);
            if untrimmed[..closing_overlap] != untrimmed[untrimmed.len() - closing_overlap..] {
                return Err("source_closing_overlap_is_not_an_exact_terminal_match".to_string());
            }
        }
    }
    if !circular && occurrences.iter().any(|copy| copy.start < 0) {
        return Err("linear_source_node_interval_precedes_path_start".to_string());
    }
    Ok(PreparedPath {
        source_spell: String::from_utf8(source).unwrap(),
        circular,
        edges: edges.to_vec(),
        occurrences,
        source_overlap_conflicts: overlap_conflicts,
    })
}

#[derive(Debug, Clone, Copy)]
struct Domain {
    low: usize,
    high: usize,
}

fn cigar_domains(source: &[u8], target: &[u8], cigar: &str) -> Result<Vec<Domain>, String> {
    let mut domains = vec![Domain { low: 0, high: 0 }; source.len() + 1];
    let (mut q, mut t, mut count) = (0usize, 0usize, 0usize);
    let mut had_digits = false;
    for op in cigar.bytes() {
        if op.is_ascii_digit() {
            count = count
                .checked_mul(10)
                .and_then(|v| v.checked_add((op - b'0') as usize))
                .ok_or("alignment_cigar_length_overflow")?;
            had_digits = true;
            continue;
        }
        if !had_digits || count == 0 {
            return Err("invalid_alignment_cigar".to_string());
        }
        match op {
            b'M' | b'=' | b'X' => {
                if count > source.len().saturating_sub(q) || count > target.len().saturating_sub(t)
                {
                    return Err("alignment_cigar_exceeds_sequences".to_string());
                }
                for _ in 0..count {
                    if (op == b'=' && source[q] != target[t])
                        || (op == b'X' && source[q] == target[t])
                    {
                        return Err(
                            "alignment_cigar_equality_operator_conflicts_with_bases".to_string()
                        );
                    }
                    q += 1;
                    t += 1;
                    domains[q] = Domain { low: t, high: t };
                }
            }
            b'I' => {
                if count > source.len().saturating_sub(q) {
                    return Err("alignment_cigar_exceeds_source".to_string());
                }
                for _ in 0..count {
                    q += 1;
                    domains[q] = Domain { low: t, high: t };
                }
            }
            b'D' => {
                if count > target.len().saturating_sub(t) {
                    return Err("alignment_cigar_exceeds_target".to_string());
                }
                t += count;
                domains[q].high = t;
            }
            _ => return Err("alignment_must_be_global_and_use_only_M_equal_X_I_D".to_string()),
        }
        count = 0;
        had_digits = false;
    }
    if had_digits || q != source.len() || t != target.len() {
        return Err("alignment_does_not_cover_both_complete_sequences".to_string());
    }
    // The frame endpoints are fixed. Terminal insertions belong to the terminal node.
    domains[0] = Domain { low: 0, high: 0 };
    domains[source.len()] = Domain {
        low: target.len(),
        high: target.len(),
    };
    Ok(domains)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AlignmentOp {
    count: usize,
    code: u8,
}

fn push_alignment_op(ops: &mut Vec<AlignmentOp>, count: usize, code: u8) {
    if count == 0 {
        return;
    }
    if let Some(last) = ops.last_mut() {
        if last.code == code {
            last.count += count;
            return;
        }
    }
    ops.push(AlignmentOp { count, code });
}

/// The caller has already validated all coordinates and equality operators.
/// Split M runs so gaps can move only across bases that really match.
fn exact_alignment_ops(source: &[u8], target: &[u8], cigar: &str) -> Vec<AlignmentOp> {
    let (mut q, mut t, mut count) = (0, 0, 0);
    let mut ops = Vec::new();
    for code in cigar.bytes() {
        if code.is_ascii_digit() {
            count = count * 10 + (code - b'0') as usize;
            continue;
        }
        match code {
            b'M' | b'=' | b'X' => {
                for _ in 0..count {
                    push_alignment_op(
                        &mut ops,
                        1,
                        if source[q] == target[t] { b'=' } else { b'X' },
                    );
                    q += 1;
                    t += 1;
                }
            }
            b'I' => {
                push_alignment_op(&mut ops, count, code);
                q += count;
            }
            b'D' => {
                push_alignment_op(&mut ops, count, code);
                t += count;
            }
            _ => unreachable!("validated alignment operation"),
        }
        count = 0;
    }
    ops
}

fn alignment_text(ops: &[AlignmentOp]) -> String {
    ops.iter()
        .map(|op| format!("{}{}", op.count, op.code as char))
        .collect()
}

const MAX_EQUIVALENT_ALIGNMENTS: usize = 256;

/// Move one gap through adjacent exact matches. These are complete equivalent
/// alignments, not independently widened cut ranges: several node boundaries must
/// remain consistent with the same actual placement of each gap.
fn equivalent_boundary_alignments(
    source: &[u8],
    target: &[u8],
    ops: &[AlignmentOp],
    boundaries: &[usize],
) -> (Vec<Vec<AlignmentOp>>, bool) {
    let mut result = Vec::new();
    let mut truncated = false;
    let (mut q, mut t) = (0usize, 0usize);
    for (index, op) in ops.iter().enumerate() {
        if matches!(op.code, b'I' | b'D') {
            let before = if index > 0 && ops[index - 1].code == b'=' {
                ops[index - 1].count
            } else {
                0
            };
            let after = if index + 1 < ops.len() && ops[index + 1].code == b'=' {
                ops[index + 1].count
            } else {
                0
            };
            let (sequence, at) = if op.code == b'I' {
                (source, q)
            } else {
                (target, t)
            };
            let mut left = 0;
            while left < before && sequence[at - left - 1] == sequence[at + op.count - left - 1] {
                left += 1;
            }
            let mut right = 0;
            while right < after && sequence[at + right] == sequence[at + op.count + right] {
                right += 1;
            }
            let mut offsets = BTreeSet::new();
            for &boundary in boundaries {
                // For I, every cut inside the deleted query interval can change.
                // For D, only crossing or landing on the cut changes its domain.
                let radius = if op.code == b'I' { op.count } else { 0 } as isize;
                let low = (boundary as isize - q as isize - radius - 1).max(-(left as isize));
                let high = (boundary as isize - q as isize + 1).min(right as isize);
                if low > high {
                    continue;
                }
                // Keep the closest placements first without allocating a range
                // proportional to a very long homopolymer or gap.
                for delta in (low.max(1)..=high).take(MAX_EQUIVALENT_ALIGNMENTS) {
                    offsets.insert(delta);
                }
                for delta in (low..=high.min(-1)).rev().take(MAX_EQUIVALENT_ALIGNMENTS) {
                    offsets.insert(delta);
                }
                if high.saturating_sub(low) as usize > 2 * MAX_EQUIVALENT_ALIGNMENTS {
                    truncated = true;
                }
            }
            let mut offsets = offsets.into_iter().collect::<Vec<_>>();
            offsets.sort_by_key(|delta| (delta.unsigned_abs(), *delta));
            for delta in offsets {
                if result.len() >= MAX_EQUIVALENT_ALIGNMENTS {
                    return (result, true);
                }
                let first = index - usize::from(before > 0);
                let last = index + 1 + usize::from(after > 0);
                let mut candidate = Vec::new();
                for part in &ops[..first] {
                    push_alignment_op(&mut candidate, part.count, part.code);
                }
                push_alignment_op(&mut candidate, (before as isize + delta) as usize, b'=');
                push_alignment_op(&mut candidate, op.count, op.code);
                push_alignment_op(&mut candidate, (after as isize - delta) as usize, b'=');
                for part in &ops[last..] {
                    push_alignment_op(&mut candidate, part.count, part.code);
                }
                result.push(candidate);
            }
        }
        if op.code != b'D' {
            q += op.count;
        }
        if op.code != b'I' {
            t += op.count;
        }
    }
    (result, truncated)
}

#[derive(Debug, Clone, Copy)]
struct Cut {
    variable: usize,
    turn: isize,
}

#[derive(Debug, Clone)]
struct ProjectedOccurrence {
    source: SourceOccurrence,
    start: Cut,
    end: Cut,
}

fn periodic_interval(
    target: &[u8],
    start: isize,
    end: isize,
    circular: bool,
) -> Result<Vec<u8>, String> {
    if end <= start {
        return Err("projected_node_would_be_empty".to_string());
    }
    if !circular {
        if start < 0 || end > target.len() as isize {
            return Err("linear_node_outside_target".to_string());
        }
        return Ok(target[start as usize..end as usize].to_vec());
    }
    let length = usize::try_from(end - start).map_err(|_| "projected_interval_length_overflow")?;
    let mut result = Vec::with_capacity(length);
    for offset in 0..length {
        let pos = (start + offset as isize).rem_euclid(target.len() as isize) as usize;
        result.push(target[pos]);
    }
    Ok(result)
}

fn cut_value(cut: Cut, assigned: &[Option<usize>], target_len: usize) -> Option<isize> {
    assigned[cut.variable].map(|v| v as isize + cut.turn * target_len as isize)
}

/// Build the portions determined by a partial boundary assignment and reject copy or
/// edge conflicts immediately. A complete successful assignment is replayed below.
fn evaluate(
    prepared: &PreparedPath,
    occurrences: &[ProjectedOccurrence],
    target: &[u8],
    assigned: &[Option<usize>],
) -> Result<Option<PathSequencePlan>, String> {
    let mut sequences: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    let mut intervals = Vec::new();
    let mut complete = true;
    for occurrence in occurrences {
        let start = cut_value(occurrence.start, assigned, target.len());
        let end = cut_value(occurrence.end, assigned, target.len());
        let interval = start.zip(end);
        intervals.push(interval);
        if let Some((start, end)) = interval {
            let sequence = oriented_bytes(
                &periodic_interval(target, start, end, prepared.circular)?,
                occurrence.source.step.reverse,
            );
            if let Some(previous) = sequences.get(&occurrence.source.step.node) {
                if *previous != sequence {
                    return Err(format!(
                        "shared_node_copies_disagree_at_projected_boundaries:{}",
                        occurrence.source.step.node
                    ));
                }
            } else {
                sequences.insert(occurrence.source.step.node.clone(), sequence);
            }
        } else {
            complete = false;
        }
    }
    let mut edge_lengths = BTreeMap::new();
    let walked_edges = occurrences
        .iter()
        .flat_map(|copy| copy.source.incoming_edges.iter().copied())
        .collect::<BTreeSet<_>>();
    for (i, occurrence) in occurrences.iter().enumerate() {
        if occurrence.source.incoming_edges.is_empty() {
            continue;
        }
        let prev = (i + occurrences.len() - 1) % occurrences.len();
        if let (Some((start, _)), Some((_, previous_end))) = (intervals[i], intervals[prev]) {
            let previous_end = previous_end - if i == 0 { target.len() as isize } else { 0 };
            let overlap = usize::try_from(previous_end - start)
                .map_err(|_| "projected_walk_has_gap_between_nodes")?;
            for edge in &occurrence.source.incoming_edges {
                if let Some(old) = edge_lengths.insert(*edge, overlap) {
                    if old != overlap {
                        return Err("same_edge_requires_different_projected_overlaps".to_string());
                    }
                }
            }
        }
    }
    for (index, edge) in prepared.edges.iter().enumerate() {
        let n = if let Some(n) = edge_lengths.get(&index) {
            *n
        } else if walked_edges.contains(&index) {
            continue;
        } else {
            edge.overlap
        };
        if let (Some(left), Some(right)) = (sequences.get(&edge.from), sequences.get(&edge.to)) {
            let left = oriented_bytes(left, edge.from_reverse);
            let right = oriented_bytes(right, edge.to_reverse);
            if n > left.len() || n > right.len() || left[left.len() - n..] != right[..n] {
                return Err(format!(
                    "edge_overlap_incompatible_with_target_sequences:{}:{}:{}:{}",
                    index, edge.from, edge.to, n
                ));
            }
        }
    }
    if !complete {
        return Ok(None);
    }
    let mut replay = Vec::new();
    let mut copies = Vec::new();
    for (i, occurrence) in occurrences.iter().enumerate() {
        let sequence = oriented_bytes(
            &sequences[&occurrence.source.step.node],
            occurrence.source.step.reverse,
        );
        let overlap = occurrence
            .source
            .incoming_edges
            .first()
            .map_or(0, |edge| edge_lengths[edge]);
        if overlap > sequence.len() {
            return Err("projected_overlap_exceeds_node_length".to_string());
        }
        replay.extend_from_slice(&sequence[overlap..]);
        let (start, end) = intervals[i].unwrap();
        copies.push(CopyPlacement {
            node: occurrence.source.step.node.clone(),
            reverse: occurrence.source.step.reverse,
            start,
            end,
        });
    }
    if replay != target {
        return Err("projected_graph_replay_does_not_equal_authoritative_target".to_string());
    }
    Ok(Some(PathSequencePlan {
        sequences: sequences
            .into_iter()
            .map(|(name, seq)| (name, String::from_utf8(seq).unwrap()))
            .collect(),
        copies,
        edge_overlaps: prepared
            .edges
            .iter()
            .enumerate()
            .map(|(index, edge)| EdgeOverlap {
                edge_index: index,
                old_overlap: edge.overlap,
                new_overlap: edge_lengths.get(&index).copied().unwrap_or(edge.overlap),
                supported_by_walk_alignment: walked_edges.contains(&index),
            })
            .collect(),
        explored_assignments: 0,
        alignment_cigar: String::new(),
        alignment_candidates_examined: 0,
        target_frame_shift: 0,
        shared_core_trims: Vec::new(),
    }))
}

/// Retain only homologous source positions whose projected target tokens agree in
/// every copy. Different flanks may be represented by the adjacent single-copy node;
/// they are never assigned by searching for an unrelated matching substring.
///
/// NOTE (repeat boundary scope, 2026-09-12): the retained homologous core can be
/// NARROWER than the repeat-node boundary inferred independently by minimap2 asm5 /
/// blastn self-alignment in the node-projection generator. This trimming only settles
/// the whole-path reconstruction boundary here; it must NOT be read as a repeat
/// upgrade/downgrade. Both the original self-alignment span (`source_length`) and the
/// retained core interval are emitted (homologous_core_trims) so a later repeat
/// upgrade/downgrade pass can reconcile the two boundary definitions explicitly.
fn project_shared_cores(
    prepared: &PreparedPath,
    target: &[u8],
    domains: &[Domain],
) -> Result<PathSequencePlan, String> {
    project_shared_cores_with_terminal_overlap(prepared, target, domains, 0)
}

fn project_shared_cores_with_terminal_overlap(
    prepared: &PreparedPath,
    target: &[u8],
    domains: &[Domain],
    terminal_overlap: usize,
) -> Result<PathSequencePlan, String> {
    let source_len = prepared.source_spell.len() as isize;
    let target_len = target.len() as isize;
    let point = |q: isize| -> Result<isize, String> {
        if terminal_overlap > 0 && q == source_len {
            return Ok(domains[q as usize].high as isize);
        }
        let (q, turn) = if prepared.circular {
            (q.rem_euclid(source_len), q.div_euclid(source_len))
        } else if (0..=source_len).contains(&q) {
            (q, 0)
        } else {
            return Err("shared_core_source_point_outside_linear_alignment".to_string());
        };
        Ok(domains[q as usize].high as isize + turn * target_len)
    };
    let mut groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    let mut intervals = Vec::new();
    for (index, copy) in prepared.occurrences.iter().enumerate() {
        groups.entry(copy.step.node.clone()).or_default().push(index);
        intervals.push((point(copy.start)?, point(copy.end)?));
    }
    let mut trims = Vec::new();
    for (node, indices) in &groups {
        if indices.len() == 1 {
            continue;
        }
        let first = &prepared.occurrences[indices[0]];
        let source_length = (first.end - first.start) as usize;
        let full = indices.iter().map(|index| {
            let copy = &prepared.occurrences[*index];
            let (start, end) = intervals[*index];
            Ok(oriented_bytes(&periodic_interval(target, start, end, prepared.circular)?, copy.step.reverse))
        }).collect::<Result<Vec<_>, String>>()?;
        if full.iter().all(|copy| *copy == full[0]) {
            continue;
        }
        let cuts = indices.iter().map(|index| {
            let copy = &prepared.occurrences[*index];
            if (copy.end - copy.start) as usize != source_length {
                return Err("shared_source_node_lengths_disagree".to_string());
            }
            (0..=source_length).map(|offset| {
                point(if copy.step.reverse { copy.end - offset as isize } else { copy.start + offset as isize })
            }).collect::<Result<Vec<_>, String>>()
        }).collect::<Result<Vec<_>, String>>()?;
        let mut run_start = 0usize;
        let mut run_bases = 0usize;
        let mut best = (0usize, 0usize, 0usize); // target bases, source start, source end
        let mut equally_best_cores = 0usize;
        for offset in 0..source_length {
            let mut reference: Option<Vec<u8>> = None;
            let mut same = true;
            for (ordinal, index) in indices.iter().enumerate() {
                let reverse = prepared.occurrences[*index].step.reverse;
                let a = cuts[ordinal][offset];
                let b = cuts[ordinal][offset + 1];
                let token = if a == b { Vec::new() } else {
                    oriented_bytes(&periodic_interval(target, a.min(b), a.max(b), prepared.circular)?, reverse)
                };
                if let Some(previous) = &reference {
                    same &= *previous == token;
                } else {
                    reference = Some(token);
                }
            }
            if same {
                run_bases += reference.as_ref().map_or(0, Vec::len);
                // Prefer the greatest retained target length, then the widest source
                // interval, then the earliest source coordinate. No free-site matching.
                let score = (run_bases, offset + 1 - run_start);
                if score > (best.0, best.2 - best.1) {
                    best = (run_bases, run_start, offset + 1);
                    equally_best_cores = 1;
                } else if score == (best.0, best.2 - best.1) {
                    equally_best_cores += 1;
                }
            } else {
                run_start = offset + 1;
                run_bases = 0;
            }
        }
        if best.0 == 0 {
            return Err(format!("shared_node_has_no_nonempty_homologous_core:{node}"));
        }
        if terminal_overlap > 0 && equally_best_cores != 1 {
            return Err(format!("circular_seam_shared_core_is_not_unique:{node}"));
        }
        for (ordinal, index) in indices.iter().enumerate() {
            let a = cuts[ordinal][best.1];
            let b = cuts[ordinal][best.2];
            intervals[*index] = (a.min(b), a.max(b));
        }
        trims.push(SharedCoreTrim {
            node: node.clone(), source_length, source_start: best.1,
            source_end: best.2, new_length: best.0,
        });
    }
    if trims.is_empty() {
        return Err("shared_core_refinement_did_not_change_a_boundary".to_string());
    }
    for index in 0..intervals.len() {
        if index == 0 && !prepared.circular {
            continue;
        }
        let previous = (index + intervals.len() - 1) % intervals.len();
        let turn = if index == 0 { target_len } else { 0 };
        let previous_end = intervals[previous].1 - turn;
        if previous_end >= intervals[index].0 {
            continue; // Existing overlap already covers the removed flank.
        }
        let left_single = groups[&prepared.occurrences[previous].step.node].len() == 1;
        let right_single = groups[&prepared.occurrences[index].step.node].len() == 1;
        if left_single {
            intervals[previous].1 = intervals[index].0 + turn;
        } else if right_single {
            intervals[index].0 = previous_end;
        } else {
            return Err("shared_core_flank_has_no_adjacent_single_copy_owner".to_string());
        }
    }
    let shift = if prepared.circular {
        intervals.last().unwrap().1 - target_len
    } else {
        if intervals[0].0 != 0 || intervals.last().unwrap().1 != target_len {
            return Err("shared_core_would_discard_a_linear_terminal_flank".to_string());
        }
        0
    };
    let frame = periodic_interval(target, shift, shift + target_len, prepared.circular)?;
    let mut assigned = Vec::new();
    let mut projected = Vec::new();
    for (copy, (start, end)) in prepared.occurrences.iter().zip(intervals) {
        let mut cut = |value: isize| {
            let value = value - shift;
            let turn = value.div_euclid(target_len);
            let variable = assigned.len();
            assigned.push(Some(value.rem_euclid(target_len) as usize));
            Cut { variable, turn }
        };
        projected.push(ProjectedOccurrence { source: copy.clone(), start: cut(start), end: cut(end) });
    }
    let mut plan = evaluate(prepared, &projected, &frame, &assigned)?
        .ok_or("shared_core_projection_is_incomplete")?;
    plan.explored_assignments = 1;
    plan.target_frame_shift = shift;
    plan.shared_core_trims = trims;
    Ok(plan)
}

/// Project every original node occurrence from the complete unfolded witness.
/// The closing overlap is derived from the endpoints, not imposed on source Ls.
fn project_circular_seam_overlap(
    prepared: &PreparedPath,
    candidate: &AlignmentCandidate,
) -> Result<PathSequencePlan, String> {
    let overlap = candidate.terminal_target_overlap;
    let target = dna(&candidate.target_in_path_frame)?;
    let source_len = prepared.source_spell.len();
    if !prepared.circular
        || overlap == 0
        || overlap >= target.len()
        || prepared.occurrences.first().map(|copy| copy.start) != Some(0)
        || prepared.occurrences.last().map(|copy| copy.end) != Some(source_len as isize)
        || prepared
            .occurrences
            .iter()
            .any(|copy| copy.start < 0 || copy.end > source_len as isize)
        || prepared.occurrences[0]
            .incoming_edges
            .iter()
            .any(|edge| prepared.edges[*edge].overlap != 0)
    {
        return Err("circular_seam_requires_complete_untrimmed_source_occurrences".to_string());
    }
    let unfolded = [target.as_slice(), &target[..overlap]].concat();
    let domains = cigar_domains(
        prepared.source_spell.as_bytes(),
        &unfolded,
        &candidate.cigar,
    )?;
    // A target-only insertion exactly on a node cut is not a unique assignment.
    for copy in &prepared.occurrences {
        for point in [copy.start, copy.end] {
            let domain = domains[point as usize];
            if domain.low != domain.high {
                return Err("circular_seam_has_ambiguous_node_boundary".to_string());
            }
        }
    }
    let frame = rotate(&candidate.target_in_path_frame, overlap);
    let mut assigned = Vec::new();
    let mut projected = Vec::new();
    for copy in &prepared.occurrences {
        let mut cut = |point: isize| {
            let value = domains[point as usize].high as isize - overlap as isize;
            let variable = assigned.len();
            assigned.push(Some(value.rem_euclid(target.len() as isize) as usize));
            Cut {
                variable,
                turn: value.div_euclid(target.len() as isize),
            }
        };
        projected.push(ProjectedOccurrence {
            source: copy.clone(),
            start: cut(copy.start),
            end: cut(copy.end),
        });
    }
    let mut plan = match evaluate(prepared, &projected, frame.as_bytes(), &assigned) {
        Ok(Some(mut plan)) => {
            plan.target_frame_shift = overlap as isize;
            plan
        }
        Err(error) if error.contains("shared_node_copies_disagree") => {
            project_shared_cores_with_terminal_overlap(prepared, &target, &domains, overlap)?
        }
        Err(error) => return Err(error),
        Ok(None) => return Err("circular_seam_projection_is_incomplete".to_string()),
    };
    plan.alignment_cigar = candidate.cigar.clone();
    plan.alignment_candidates_examined = 1;
    plan.explored_assignments = 1;
    Ok(plan)
}

#[derive(Debug)]
struct BoundaryFailure {
    reason: String,
    explored: usize,
}

/// Construct a sequence plan from one complete alignment. At a source boundary with
/// target-only insertion, every position inside that insertion is an evidenced cut.
/// Shared nodes constrain those cuts jointly. The default assigns an insertion to the
/// preceding interval; alternatives are visited only when the joint constraints fail.
/// If shared copies still disagree, gaps may move across adjacent exact matches.
/// Every such alternative remains one complete alignment with unchanged edit cost;
/// its node cuts are solved together and the whole graph is checked again.
///
/// `search_limit` is a computational budget, not an acceptance threshold. Exhaustion or
/// failure means this alignment did not establish a solution; it is not proof that no
/// topology-preserving solution exists under another equivalent alignment.
pub fn project_path_sequence(
    prepared: &PreparedPath,
    target_in_path_frame: &str,
    global_cigar: &str,
    search_limit: usize,
) -> Result<PathSequencePlan, String> {
    if search_limit == 0 {
        return Err("boundary_search_limit_must_be_positive".to_string());
    }
    let target = dna(target_in_path_frame)?;
    let full_domains = cigar_domains(prepared.source_spell.as_bytes(), &target, global_cigar)?;
    let first = project_boundary_domains(prepared, &target, &full_domains, search_limit);
    let failure = match first {
        Ok(mut plan) => {
            plan.alignment_cigar = global_cigar.to_string();
            plan.alignment_candidates_examined = 1;
            return Ok(plan);
        }
        Err(failure) => failure,
    };
    // A disagreement under one alignment is not evidence that the repeat copies
    // differ. A gap next to a shared boundary may have an equivalent placement.
    if !failure.reason.contains("shared_node_copies_disagree") {
        return Err(failure.reason);
    }
    let source = prepared.source_spell.as_bytes();
    let mut boundaries = BTreeSet::new();
    for copy in &prepared.occurrences {
        for point in [copy.start, copy.end] {
            boundaries.insert(if prepared.circular {
                point.rem_euclid(source.len() as isize) as usize
            } else {
                point as usize
            });
        }
    }
    let boundaries = boundaries.into_iter().collect::<Vec<_>>();
    let original = exact_alignment_ops(source, &target, global_cigar);
    let mut seen = BTreeSet::from([alignment_text(&original)]);
    let mut pending = VecDeque::from([original]);
    let mut total_explored = failure.explored;
    let mut examined = 1;
    let mut truncated = false;
    let mut last_conflict = failure.reason;
    // Keep full original node boundaries when an equivalent alignment can satisfy
    // them. A homologous-core candidate is a distinct, explicitly recorded boundary
    // refinement and still must pass the same all-edge and full-sequence checks.
    let mut core_result = if total_explored < search_limit {
        Some(project_shared_cores(prepared, &target, &full_domains))
    } else {
        None
    };
    let accept_core = |plan: &mut PathSequencePlan, explored, examined| {
        plan.explored_assignments += explored;
        plan.alignment_candidates_examined = examined;
        plan.alignment_cigar = global_cigar.to_string();
    };
    while let Some(ops) = pending.pop_front() {
        let (candidates, cut_short) =
            equivalent_boundary_alignments(source, &target, &ops, &boundaries);
        truncated |= cut_short;
        for candidate in candidates {
            let cigar = alignment_text(&candidate);
            if !seen.insert(cigar.clone()) {
                continue;
            }
            if examined >= MAX_EQUIVALENT_ALIGNMENTS || total_explored >= search_limit {
                if let Some(Ok(mut plan)) = core_result.take() {
                    accept_core(&mut plan, total_explored, examined);
                    return Ok(plan);
                }
                return Err(format!(
                    "equivalent_alignment_search_resource_limit:alignments={examined}:assignments={total_explored}:{last_conflict}"
                ));
            }
            // Validate every actual witness; moving a gap can never turn a
            // mismatch into a match or hide an uncovered terminal sequence.
            let domains = cigar_domains(source, &target, &cigar)?;
            examined += 1;
            match project_boundary_domains(
                prepared,
                &target,
                &domains,
                search_limit - total_explored,
            ) {
                Ok(mut plan) => {
                    plan.explored_assignments += total_explored;
                    plan.alignment_candidates_examined = examined;
                    plan.alignment_cigar = cigar;
                    return Ok(plan);
                }
                Err(error) => {
                    total_explored += error.explored;
                    last_conflict = error.reason;
                }
            }
            pending.push_back(candidate);
        }
    }
    match core_result {
        Some(Ok(mut plan)) => {
            accept_core(&mut plan, total_explored, examined);
            return Ok(plan);
        }
        Some(Err(reason)) => last_conflict.push_str(&format!(";shared_core_refinement:{reason}")),
        None => {}
    }
    let reason = if truncated {
        "equivalent_alignment_search_resource_limit"
    } else {
        "equivalent_alignment_boundary_constraints_unresolved"
    };
    Err(format!(
        "{reason}:alignments={examined}:assignments={total_explored}:{last_conflict}"
    ))
}

fn project_boundary_domains(
    prepared: &PreparedPath,
    target: &[u8],
    full_domains: &[Domain],
    search_limit: usize,
) -> Result<PathSequencePlan, BoundaryFailure> {
    let source_len = prepared.source_spell.len() as isize;
    let mut variable_sources = BTreeSet::new();
    for copy in &prepared.occurrences {
        for point in [copy.start, copy.end] {
            if !prepared.circular && (point < 0 || point > source_len) {
                return Err(BoundaryFailure {
                    reason: "source_node_interval_outside_linear_alignment".to_string(),
                    explored: 0,
                });
            }
            variable_sources.insert(if prepared.circular {
                point.rem_euclid(source_len)
            } else {
                point
            } as usize);
        }
    }
    let variable_sources = variable_sources.into_iter().collect::<Vec<_>>();
    let variable_index = variable_sources
        .iter()
        .enumerate()
        .map(|(i, q)| (*q, i))
        .collect::<BTreeMap<_, _>>();
    let domains = variable_sources
        .iter()
        .map(|q| full_domains[*q])
        .collect::<Vec<_>>();
    let cut = |point: isize| -> Cut {
        let (q, turn) = if prepared.circular {
            (point.rem_euclid(source_len), point.div_euclid(source_len))
        } else {
            (point, 0)
        };
        Cut {
            variable: variable_index[&(q as usize)],
            turn,
        }
    };
    let occurrences = prepared
        .occurrences
        .iter()
        .map(|copy| ProjectedOccurrence {
            source: copy.clone(),
            start: cut(copy.start),
            end: cut(copy.end),
        })
        .collect::<Vec<_>>();
    let mut assigned = domains
        .iter()
        .map(|domain| (domain.low == domain.high).then_some(domain.low))
        .collect::<Vec<_>>();
    let mut multiplicity = BTreeMap::new();
    for copy in &occurrences {
        *multiplicity
            .entry(copy.source.step.node.clone())
            .or_insert(0usize) += 1;
    }
    let mut priority = vec![0usize; domains.len()];
    for copy in &occurrences {
        let weight = multiplicity[&copy.source.step.node];
        priority[copy.start.variable] += weight * weight;
        priority[copy.end.variable] += weight * weight;
    }
    let mut variables = (0..domains.len())
        .filter(|i| assigned[*i].is_none())
        .collect::<Vec<_>>();
    variables.sort_by_key(|i| {
        (
            std::cmp::Reverse(priority[*i]),
            domains[*i].high - domains[*i].low,
            *i,
        )
    });
    struct Search<'a> {
        prepared: &'a PreparedPath,
        occurrences: &'a [ProjectedOccurrence],
        target: &'a [u8],
        domains: &'a [Domain],
        variables: &'a [usize],
        limit: usize,
        explored: usize,
        last_conflict: String,
    }
    fn visit(
        search: &mut Search<'_>,
        assigned: &mut [Option<usize>],
        depth: usize,
    ) -> Option<PathSequencePlan> {
        if search.explored >= search.limit {
            return None;
        }
        search.explored += 1;
        match evaluate(search.prepared, search.occurrences, search.target, assigned) {
            Ok(Some(mut plan)) => {
                plan.explored_assignments = search.explored;
                return Some(plan);
            }
            Err(conflict) => {
                search.last_conflict = conflict;
                return None;
            }
            Ok(None) => {}
        }
        let variable = *search.variables.get(depth)?;
        let domain = search.domains[variable];
        for value in (domain.low..=domain.high).rev() {
            assigned[variable] = Some(value);
            if let Some(plan) = visit(search, assigned, depth + 1) {
                return Some(plan);
            }
            if search.explored >= search.limit {
                break;
            }
        }
        assigned[variable] = None;
        None
    }
    let mut search = Search {
        prepared,
        occurrences: &occurrences,
        target,
        domains: &domains,
        variables: &variables,
        limit: search_limit,
        explored: 0,
        last_conflict: "no_complete_boundary_assignment".to_string(),
    };
    if let Some(plan) = visit(&mut search, &mut assigned, 0) {
        return Ok(plan);
    }
    let reason = if search.explored >= search.limit {
        "boundary_search_resource_limit"
    } else {
        "alignment_boundary_constraints_unresolved"
    };
    Err(BoundaryFailure {
        reason: format!(
            "{}:assignments={}:{}",
            reason, search.explored, search.last_conflict
        ),
        explored: search.explored,
    })
}

#[cfg(test)]
mod tests_path_sequence {
    use super::*;

    fn step(node: &str, reverse: bool) -> PathStep {
        PathStep {
            node: node.into(),
            reverse,
        }
    }
    fn edge(a: &str, ar: bool, b: &str, br: bool, overlap: usize) -> PathEdge {
        PathEdge {
            from: a.into(),
            from_reverse: ar,
            to: b.into(),
            to_reverse: br,
            overlap,
        }
    }
    fn nodes(values: &[(&str, &str)]) -> BTreeMap<String, String> {
        values
            .iter()
            .map(|(name, seq)| (name.to_string(), seq.to_string()))
            .collect()
    }

    #[test]
    fn target_insertion_at_zero_overlap_is_assigned_without_a_gap() {
        let p = prepare_path(
            &nodes(&[("a", "AA"), ("b", "CC")]),
            &[edge("a", false, "b", false, 0)],
            &[step("a", false), step("b", false)],
            false,
        )
        .unwrap();
        let plan = project_path_sequence(&p, "AAGGCC", "2=2D2=", 100).unwrap();
        assert_eq!(plan.sequences["a"], "AAGG");
        assert_eq!(plan.sequences["b"], "CC");
        assert_eq!(plan.edge_overlaps[0].new_overlap, 0);
    }

    #[test]
    fn an_insertion_inside_an_overlap_updates_both_node_copies_and_edge_length() {
        let p = prepare_path(
            &nodes(&[("a", "AACG"), ("b", "CGTT")]),
            &[edge("a", false, "b", false, 2)],
            &[step("a", false), step("b", false)],
            false,
        )
        .unwrap();
        assert_eq!(p.source_spell, "AACGTT");
        let plan = project_path_sequence(&p, "AACAGTT", "3=1D3=", 100).unwrap();
        assert_eq!(plan.sequences["a"], "AACAG");
        assert_eq!(plan.sequences["b"], "CAGTT");
        assert_eq!(plan.edge_overlaps[0].new_overlap, 3);
    }

    #[test]
    fn conflicting_source_overlap_cannot_be_repaired_by_discarding_a_node_prefix() {
        let error = prepare_path(
            &nodes(&[("a", "AACC"), ("b", "GGTT")]),
            &[edge("a", false, "b", false, 2)],
            &[step("a", false), step("b", false)],
            false,
        )
        .unwrap_err();
        assert!(error.contains("source_edge_overlap_sequence_conflict:a:b:2"));
    }

    #[test]
    fn reverse_oriented_shared_copies_must_be_the_same_forward_sequence() {
        let ns = nodes(&[("r", "AC"), ("b", "GG")]);
        let es = [
            edge("r", false, "b", false, 0),
            edge("b", false, "r", true, 0),
        ];
        let p = prepare_path(
            &ns,
            &es,
            &[step("r", false), step("b", false), step("r", true)],
            false,
        )
        .unwrap();
        let plan = project_path_sequence(&p, "ATGGAT", "1=1X2=1X1=", 100).unwrap();
        assert_eq!(plan.sequences["r"], "AT");
    }

    #[test]
    fn shared_copy_conflict_rejects_instead_of_splitting_the_node() {
        let ns = nodes(&[("r", "AC"), ("b", "GG")]);
        let p = prepare_path(
            &ns,
            &[
                edge("r", false, "b", false, 0),
                edge("b", false, "r", false, 0),
            ],
            &[step("r", false), step("b", false), step("r", false)],
            false,
        )
        .unwrap();
        let error = project_path_sequence(&p, "ACGGAT", "5=1X", 100).unwrap_err();
        assert!(error.contains("shared_node_copies_disagree"));
    }

    #[test]
    fn ambiguous_insertions_are_assigned_jointly_across_repeat_copies() {
        let ns = nodes(&[("r", "AC"), ("b", "GG")]);
        let p = prepare_path(
            &ns,
            &[
                edge("r", false, "b", false, 0),
                edge("b", false, "r", false, 0),
            ],
            &[step("r", false), step("b", false), step("r", false)],
            false,
        )
        .unwrap();
        // The insertion after the first r belongs to b, since the other copy ends
        // at the fixed target endpoint and has no insertion.
        let plan = project_path_sequence(&p, "ACTGGAC", "2=1D4=", 100).unwrap();
        assert_eq!(plan.sequences["r"], "AC");
        assert_eq!(plan.sequences["b"], "TGG");
        assert!(plan.explored_assignments > 1);
    }

    #[test]
    fn a_deletion_at_a_repeat_boundary_can_belong_to_the_single_copy_node() {
        let ns = nodes(&[("a", "GC"), ("r", "TA"), ("b", "ATC")]);
        let es = [
            edge("a", false, "r", false, 0),
            edge("r", false, "b", false, 0),
            edge("b", false, "r", false, 0),
        ];
        let p = prepare_path(
            &ns,
            &es,
            &[
                step("a", false),
                step("r", false),
                step("b", false),
                step("r", false),
            ],
            false,
        )
        .unwrap();
        assert_eq!(p.source_spell, "GCTAATCTA");
        let plan = project_path_sequence(&p, "GCTATCTA", "3=1I5=", 100).unwrap();
        assert_eq!(plan.sequences["r"], "TA");
        assert_eq!(plan.sequences["b"], "TC");
        assert_eq!(plan.alignment_cigar, "4=1I4=");
        assert!(plan.alignment_candidates_examined > 1);
        assert!(plan.edge_overlaps.iter().all(|edge| edge.new_overlap == 0));
    }

    #[test]
    fn equivalent_indel_placement_respects_reverse_oriented_repeat_copies() {
        let ns = nodes(&[("a", "GC"), ("r", "AC"), ("b", "TGA")]);
        let es = [
            edge("a", false, "r", true, 0),
            edge("r", true, "b", false, 0),
            edge("b", false, "r", false, 0),
        ];
        let p = prepare_path(
            &ns,
            &es,
            &[
                step("a", false),
                step("r", true),
                step("b", false),
                step("r", false),
            ],
            false,
        )
        .unwrap();
        let plan = project_path_sequence(&p, "GCGTGAAC", "3=1I5=", 100).unwrap();
        assert_eq!(plan.sequences["r"], "AC");
        assert_eq!(plan.sequences["b"], "GA");
        assert!(plan.copies[1].reverse);
    }

    #[test]
    fn both_ends_of_a_single_copy_node_can_shift_without_changing_the_repeat() {
        let ns = nodes(&[("a", "GC"), ("r", "TA"), ("b", "AAG")]);
        let es = [
            edge("a", false, "r", false, 0),
            edge("r", false, "b", false, 0),
            edge("b", false, "r", false, 0),
            edge("r", false, "a", false, 0),
        ];
        let p = prepare_path(
            &ns,
            &es,
            &[
                step("a", false),
                step("r", false),
                step("b", false),
                step("r", false),
            ],
            true,
        )
        .unwrap();
        let plan = project_path_sequence(&p, "GCTAAGTTA", "3=1I4=1D1=", 100).unwrap();
        assert_eq!(plan.sequences["r"], "TA");
        assert_eq!(plan.sequences["b"], "AGT");
        assert_eq!(plan.sequences["a"], "GC");
        assert_eq!(plan.copies.len(), 4);
        assert!(plan.alignment_candidates_examined > 2);
    }

    #[test]
    fn equivalent_indel_search_still_checks_edges_outside_the_walk() {
        let ns = nodes(&[("a", "GC"), ("r", "TA"), ("b", "ATC")]);
        let es = [
            edge("a", false, "r", false, 0),
            edge("r", false, "b", false, 0),
            edge("b", false, "r", false, 0),
            edge("b", false, "a", false, 1),
        ];
        let p = prepare_path(
            &ns,
            &es,
            &[
                step("a", false),
                step("r", false),
                step("b", false),
                step("r", false),
            ],
            false,
        )
        .unwrap();
        assert!(project_path_sequence(&p, "GCTATCTA", "3=1I5=", 100).is_err());
    }

    #[test]
    fn equivalent_indel_search_respects_the_total_assignment_budget() {
        let ns = nodes(&[("a", "GC"), ("r", "TA"), ("b", "ATC")]);
        let p = prepare_path(
            &ns,
            &[
                edge("a", false, "r", false, 0),
                edge("r", false, "b", false, 0),
                edge("b", false, "r", false, 0),
            ],
            &[
                step("a", false),
                step("r", false),
                step("b", false),
                step("r", false),
            ],
            false,
        )
        .unwrap();
        let error = project_path_sequence(&p, "GCTATCTA", "3=1I5=", 1).unwrap_err();
        assert!(error.contains("equivalent_alignment_search_resource_limit"));
    }

    #[test]
    fn unused_edges_are_checked_and_never_adjusted_without_walk_evidence() {
        let ns = nodes(&[("a", "AAC"), ("b", "CGG")]);
        let es = [
            edge("a", false, "b", false, 1),
            edge("b", false, "a", false, 1),
        ];
        let p = prepare_path(&ns, &es, &[step("a", false), step("b", false)], false).unwrap();
        let error = project_path_sequence(&p, "AACGG", "5=", 100).unwrap_err();
        assert!(error.contains("edge_overlap_incompatible"));
    }

    #[test]
    fn a_closed_path_keeps_terminal_overlap_once() {
        let ns = nodes(&[("a", "ACGTAC")]);
        let p = prepare_path(
            &ns,
            &[edge("a", false, "a", false, 2)],
            &[step("a", false)],
            true,
        )
        .unwrap();
        assert_eq!(p.source_spell, "GTAC");
        let plan = project_path_sequence(&p, "GTAC", "4=", 100).unwrap();
        assert_eq!(plan.sequences["a"], "ACGTAC");
        assert_eq!(plan.copies[0].start, -2);
        assert_eq!(plan.edge_overlaps[0].new_overlap, 2);
    }

    #[test]
    fn closing_overlap_changes_are_projected_across_the_frame_origin() {
        let ns = nodes(&[("a", "ACGTAC")]);
        let p = prepare_path(
            &ns,
            &[edge("a", false, "a", false, 2)],
            &[step("a", false)],
            true,
        )
        .unwrap();
        let plan = project_path_sequence(&p, "GTAGC", "3=1D1=", 100).unwrap();
        assert_eq!(plan.sequences["a"], "AGCGTAGC");
        assert_eq!(plan.edge_overlaps[0].new_overlap, 3);
        assert_eq!(plan.copies[0].start, -3);
    }

    #[test]
    fn a_nonmatching_original_closure_cannot_be_silently_trimmed() {
        let ns = nodes(&[("a", "ACGTAG")]);
        let error = prepare_path(
            &ns,
            &[edge("a", false, "a", false, 2)],
            &[step("a", false)],
            true,
        )
        .unwrap_err();
        assert!(error.contains("source_closing_overlap_is_not_an_exact_terminal_match"));
    }

    #[test]
    fn unaligned_ends_and_wrong_direction_do_not_pass() {
        let ns = nodes(&[("a", "AAC"), ("b", "CGG")]);
        assert!(prepare_path(
            &ns,
            &[edge("a", false, "b", true, 1)],
            &[step("a", false), step("b", false)],
            false
        )
        .is_err());
        let p = prepare_path(&nodes(&[("a", "AAC")]), &[], &[step("a", false)], false).unwrap();
        assert!(project_path_sequence(&p, "AAC", "2=", 100).is_err());
    }
}



// =====================================================================
// Inlined module: rebuild_recorded_path.rs (whole-path generator, generator 1)
// =====================================================================
// Integrate the checked topology and fixed-walk sequence constraints into rebuild.


type BuildResult = (
    SequencePlan,
    BTreeMap<String, String>,
    Vec<(&'static str, String)>,
);

/// A representation preference, not a claim that an assignment is biologically unique.
/// Only complete plans already accepted by the fixed-walk sequence constraints are ranked.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct AssignmentPreference {
    changed_overlaps: usize,
    overlap_delta: u128,
    changed_nodes: usize,
    node_length_delta: u128,
    changed_sequence_span: u128,
}

fn assignment_preference(
    plan: &PathSequencePlan,
    source: &BTreeMap<String, String>,
) -> AssignmentPreference {
    let mut score = AssignmentPreference {
        changed_overlaps: 0,
        overlap_delta: 0,
        changed_nodes: 0,
        node_length_delta: 0,
        changed_sequence_span: 0,
    };
    for edge in &plan.edge_overlaps {
        if edge.old_overlap != edge.new_overlap {
            score.changed_overlaps += 1;
            score.overlap_delta += edge.old_overlap.abs_diff(edge.new_overlap) as u128;
        }
    }
    for (node, new) in &plan.sequences {
        let old = &source[node];
        if old == new {
            continue;
        }
        score.changed_nodes += 1;
        score.node_length_delta += old.len().abs_diff(new.len()) as u128;
        // This is the span outside common terminal sequence, not edit distance.
        // It favors retaining the original node origin without an expensive new alignment.
        let prefix = old
            .bytes()
            .zip(new.bytes())
            .take_while(|(a, b)| a == b)
            .count();
        let suffix = old.as_bytes()[prefix..]
            .iter()
            .rev()
            .zip(new.as_bytes()[prefix..].iter().rev())
            .take_while(|(a, b)| a == b)
            .count();
        score.changed_sequence_span += (old.len().max(new.len()) - prefix - suffix) as u128;
    }
    score
}

type AcceptedAssignment = (PathSequencePlan, usize, String, String);

fn compare_assignments(
    left: &AcceptedAssignment,
    right: &AcceptedAssignment,
    source: &BTreeMap<String, String>,
) -> std::cmp::Ordering {
    assignment_preference(&left.0, source)
        .cmp(&assignment_preference(&right.0, source))
        // Exact content is a stable tie-break, never the order minimap2 returned hits.
        .then_with(|| left.0.sequences.cmp(&right.0.sequences))
        .then_with(|| {
            left.0
                .edge_overlaps
                .iter()
                .map(|e| (e.edge_index, e.new_overlap))
                .cmp(
                    right
                        .0
                        .edge_overlaps
                        .iter()
                        .map(|e| (e.edge_index, e.new_overlap)),
                )
        })
        .then_with(|| left.1.cmp(&right.1))
        .then_with(|| {
            left.0
                .copies
                .iter()
                .map(|c| (&c.node, c.reverse, c.start, c.end))
                .cmp(
                    right
                        .0
                        .copies
                        .iter()
                        .map(|c| (&c.node, c.reverse, c.start, c.end)),
                )
        })
        .then_with(|| left.3.cmp(&right.3))
}

fn build_plan(
    graph: &mut Gfa,
    target: &FastaRecord,
    contract: &ResolveContract,
    minimap2: &Path,
    threads: usize,
    paths: &OutputPaths,
) -> Result<BuildResult, OrgraftError> {
    let invalid = |message: String| OrgraftError::InvalidArgument(message);
    if target.sequence.is_empty() {
        return Err(invalid("authoritative FASTA is empty".to_string()));
    }
    let edges = contract
        .topology_links
        .iter()
        .map(|edge| {
            Ok(PathEdge {
                from: edge.from.clone(),
                from_reverse: edge.from_orientation == '-',
                to: edge.to.clone(),
                to_reverse: edge.to_orientation == '-',
                overlap: edge.overlap_bases()?,
            })
        })
        .collect::<Result<Vec<_>, OrgraftError>>()?;
    let walk = contract
        .walk
        .iter()
        .map(|step| PathStep {
            node: step.source_node.clone(),
            reverse: step.orientation == '-',
        })
        .collect::<Vec<_>>();
    // The whole-path generator tolerates documented source overlap end-sequence errors
    // (polish corrected them) and records each conflict; finish_plan still requires every
    // output edge to be exactly compatible and the full replay to equal the target.
    let prepared = prepare_path_corrected(
        &contract.source_sequences,
        &edges,
        &walk,
        contract.circular,
    )
    .map_err(invalid)?;
    let source_overlap_conflicts = prepared.source_overlap_conflicts.clone();
    // Alignment evidence is kept with the run, including failed attempts.
    fs::create_dir_all(&paths.debug_dir)?;
    let candidates = align_path(
        &prepared.source_spell,
        &target.sequence,
        contract.circular,
        contract.alignment_reverse_complemented,
        minimap2,
        threads,
        &paths.debug_dir,
    )
    .map_err(invalid)?;
    let mut accepted: Option<AcceptedAssignment> = None;
    let mut multiple_valid_assignments = false;
    let mut failures = Vec::new();
    let mut accepted_count = 0usize;
    for candidate in candidates {
        let projected = if candidate.terminal_target_overlap > 0 {
            project_circular_seam_overlap(&prepared, &candidate)
        } else {
            project_path_sequence(
                &prepared,
                &candidate.target_in_path_frame,
                &candidate.cigar,
                10_000,
            )
        };
        match projected {
            Ok(plan) => {
                let shift = plan.target_frame_shift.rem_euclid(target.sequence.len() as isize) as usize;
                let frame = format!("{}{}", &candidate.target_in_path_frame[shift..], &candidate.target_in_path_frame[..shift]);
                let proposed = (
                    plan,
                    (candidate.frame_offset + shift) % target.sequence.len(),
                    frame,
                    candidate.method,
                );
                if let Some(previous) = &accepted {
                    multiple_valid_assignments |= previous.0.sequences != proposed.0.sequences
                        || previous.0.edge_overlaps != proposed.0.edge_overlaps;
                    if compare_assignments(&proposed, previous, &contract.source_sequences).is_lt()
                    {
                        accepted = Some(proposed);
                    }
                } else {
                    accepted = Some(proposed);
                }
                accepted_count += 1;
            }
            Err(error) => failures.push(error),
        }
    }
    let assignment = accepted.ok_or_else(|| {
        invalid(format!(
            "no sequence assignment was established for the recorded path; {}",
            failures.join("; ")
        ))
    })?;
    finish_plan(
        graph,
        target,
        contract,
        paths,
        assignment,
        CompletionEvidence {
            source_path_length: Some(prepared.source_spell.len()),
            accepted_alignments: accepted_count,
            multiple_valid_assignments,
            method: "recorded_resolve_path_alignment".to_string(),
            source_overlap_conflicts,
        },
    )
}

#[derive(Debug)]
struct CompletionEvidence {
    source_path_length: Option<usize>,
    accepted_alignments: usize,
    multiple_valid_assignments: bool,
    method: String,
    source_overlap_conflicts: Vec<SourceOverlapConflict>,
}

fn oriented_endpoints(
    from: &str,
    from_reverse: bool,
    to: &str,
    to_reverse: bool,
) -> ((String, bool), (String, bool)) {
    let a = (from.to_string(), from_reverse);
    let b = (to.to_string(), !to_reverse);
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

/// Keep edge indices tied to the approved contract; overlap lengths may be
/// projected, but no caller may supply a different directed connection.
fn contract_edges(
    graph: &Gfa,
    contract: &ResolveContract,
) -> Result<Vec<PathEdge>, OrgraftError> {
    let invalid = |message: &str| OrgraftError::InvalidArgument(message.to_string());
    let names = graph.order.iter().cloned().collect::<BTreeSet<_>>();
    if names.len() != graph.order.len()
        || names != graph.segments.keys().cloned().collect()
        || names != contract.source_sequences.keys().cloned().collect()
    {
        return Err(invalid(
            "candidate graph node identities differ from the approved contract",
        ));
    }
    if graph.links.len() != contract.topology_links.len() {
        return Err(invalid(
            "candidate graph link count differs from the approved contract",
        ));
    }
    let mut edges = Vec::new();
    for (actual, approved) in graph.links.iter().zip(&contract.topology_links) {
        if !matches!(actual.from_orient, '+' | '-') || !matches!(actual.to_orient, '+' | '-') {
            return Err(invalid("candidate graph has an invalid edge orientation"));
        }
        let actual_key = oriented_endpoints(
            &actual.from_name,
            actual.from_orient == '-',
            &actual.to_name,
            actual.to_orient == '-',
        );
        let approved_key = oriented_endpoints(
            &approved.from,
            approved.from_orientation == '-',
            &approved.to,
            approved.to_orientation == '-',
        );
        if actual_key != approved_key {
            return Err(invalid(
                "candidate graph directed connections differ from the approved contract",
            ));
        }
        let checked = EvidenceLink {
            from: actual.from_name.clone(),
            from_orientation: actual.from_orient,
            to: actual.to_name.clone(),
            to_orientation: actual.to_orient,
            overlap: actual.overlap.clone(),
            tags: Vec::new(),
        };
        edges.push(PathEdge {
            from: actual.from_name.clone(),
            from_reverse: actual.from_orient == '-',
            to: actual.to_name.clone(),
            to_reverse: actual.to_orient == '-',
            overlap: checked.overlap_bases()?,
        });
    }
    Ok(edges)
}

fn all_edge_overlaps(
    sequences: &BTreeMap<String, String>,
    edges: &[PathEdge],
) -> Result<(), OrgraftError> {
    for edge in edges {
        let left = sequences.get(&edge.from).ok_or_else(|| {
            OrgraftError::InvalidArgument("candidate edge source node is missing".to_string())
        })?;
        let right = sequences.get(&edge.to).ok_or_else(|| {
            OrgraftError::InvalidArgument("candidate edge target node is missing".to_string())
        })?;
        let left = if edge.from_reverse {
            reverse_complement(left)
        } else {
            left.clone()
        };
        let right = if edge.to_reverse {
            reverse_complement(right)
        } else {
            right.clone()
        };
        let n = edge.overlap;
        if n > left.len()
            || n > right.len()
            || left.as_bytes()[left.len() - n..] != right.as_bytes()[..n]
        {
            return Err(OrgraftError::InvalidArgument(format!(
                "candidate edge overlap conflicts with node sequences:{}:{}:{}",
                edge.from, edge.to, n
            )));
        }
    }
    Ok(())
}

/// Accept freshly inferred node sequences only after the same complete graph and
/// authoritative FASTA checks used for whole-path alignment plans. This function
/// reads no previous model and writes no alignment evidence it did not compute.
fn accept_node_sequences(
    graph: &mut Gfa,
    target: &FastaRecord,
    contract: &ResolveContract,
    sequences: BTreeMap<String, String>,
    paths: &OutputPaths,
    method: &str,
) -> Result<BuildResult, OrgraftError> {
    let invalid = |message: &str| OrgraftError::InvalidArgument(message.to_string());
    if target.sequence.is_empty() {
        return Err(invalid("authoritative FASTA is empty"));
    }
    if sequences.keys().cloned().collect::<BTreeSet<_>>()
        != contract.source_sequences.keys().cloned().collect()
    {
        return Err(invalid(
            "candidate node sequences change node identity or shared repeat representation",
        ));
    }
    let sequences = sequences
        .into_iter()
        .map(|(node, sequence)| (node, sequence.to_ascii_uppercase()))
        .collect::<BTreeMap<_, _>>();
    let edges = contract_edges(graph, contract)?;
    let walk = contract
        .walk
        .iter()
        .map(|step| PathStep {
            node: step.source_node.clone(),
            reverse: step.orientation == '-',
        })
        .collect::<Vec<_>>();
    all_edge_overlaps(&sequences, &edges)?;
    let prepared =
        prepare_path(&sequences, &edges, &walk, contract.circular)
            .map_err(OrgraftError::InvalidArgument)?;
    let oriented_target = if contract.alignment_reverse_complemented {
        reverse_complement(&target.sequence.to_ascii_uppercase())
    } else {
        target.sequence.to_ascii_uppercase()
    };
    if prepared.source_spell.len() != oriented_target.len() {
        return Err(invalid(
            "candidate graph replay length differs from the authoritative FASTA",
        ));
    }
    let offset = if contract.circular {
        format!("{oriented_target}{oriented_target}").find(&prepared.source_spell).filter(|offset| *offset < oriented_target.len())
    } else if prepared.source_spell == oriented_target { Some(0) } else { None }
        .ok_or_else(|| invalid("candidate graph replay does not equal the authoritative FASTA in the allowed direction and origin"))?;
    let source_path_length = prepare_path(
        &contract.source_sequences,
        &edges,
        &walk,
        contract.circular,
    )
    .ok()
    .map(|source| source.source_spell.len());
    let plan = PathSequencePlan {
        sequences,
        copies: prepared.source_copies(),
        edge_overlaps: edges
            .iter()
            .enumerate()
            .map(|(index, edge)| {
                Ok(EdgeOverlap {
                    edge_index: index,
                    old_overlap: contract.topology_links[index].overlap_bases()?,
                    new_overlap: edge.overlap,
                    supported_by_walk_alignment: false,
                })
            })
            .collect::<Result<Vec<_>, OrgraftError>>()?,
        explored_assignments: 0,
        alignment_cigar: ".".to_string(),
        alignment_candidates_examined: 0,
        target_frame_shift: 0,
        shared_core_trims: Vec::new(),
    };
    finish_plan(
        graph,
        target,
        contract,
        paths,
        (plan, offset, prepared.source_spell, ".".to_string()),
        CompletionEvidence {
            source_path_length,
            accepted_alignments: 0,
            multiple_valid_assignments: false,
            method: method.to_string(),
            source_overlap_conflicts: Vec::new(),
        },
    )
}

fn finish_plan(
    graph: &mut Gfa,
    target: &FastaRecord,
    contract: &ResolveContract,
    paths: &OutputPaths,
    assignment: AcceptedAssignment,
    evidence: CompletionEvidence,
) -> Result<BuildResult, OrgraftError> {
    let invalid = |message: String| OrgraftError::InvalidArgument(message);
    let (plan, frame_offset, frame, alignment_method) = assignment;
    let edges = contract_edges(graph, contract)?;
    let walk = contract
        .walk
        .iter()
        .map(|step| PathStep {
            node: step.source_node.clone(),
            reverse: step.orientation == '-',
        })
        .collect::<Vec<_>>();
    let selected_preference = assignment_preference(&plan, &contract.source_sequences);
    let mut updated_edges = edges.clone();
    let mut updated_links = graph.links.clone();
    if plan.edge_overlaps.len() != edges.len()
        || plan
            .edge_overlaps
            .iter()
            .enumerate()
            .any(|(index, edge)| index != edge.edge_index)
    {
        return Err(invalid(
            "candidate overlap plan does not cover every approved graph edge".to_string(),
        ));
    }
    for overlap in &plan.edge_overlaps {
        let link = &mut updated_links[overlap.edge_index];
        if overlap.old_overlap != overlap.new_overlap {
            link.tags = replace_tag(&link.tags, format!("PO:Z:{}M", overlap.old_overlap));
            link.overlap = format!("{}M", overlap.new_overlap);
        }
        updated_edges[overlap.edge_index].overlap = overlap.new_overlap;
    }
    all_edge_overlaps(&plan.sequences, &updated_edges)?;
    // Reconstruct again from the final node sequences and updated links. This cannot
    // pass merely because the output FASTA is a copy of the authoritative input.
    let replay = prepare_path(
        &plan.sequences,
        &updated_edges,
        &walk,
        contract.circular,
    )
    .map_err(invalid)?;
    if replay.source_spell != frame {
        return Err(invalid(
            "final graph replay differs from the authoritative target frame".to_string(),
        ));
    }
    let oriented_target = if contract.alignment_reverse_complemented {
        reverse_complement(&target.sequence.to_ascii_uppercase())
    } else {
        target.sequence.to_ascii_uppercase()
    };
    if oriented_target.is_empty()
        || frame.len() != oriented_target.len()
        || frame_offset >= oriented_target.len()
        || (!contract.circular && frame_offset != 0)
        || frame
            != format!(
                "{}{}",
                &oriented_target[frame_offset..],
                &oriented_target[..frame_offset]
            )
    {
        return Err(invalid("accepted graph frame does not reproduce the authoritative FASTA under its allowed origin".to_string()));
    }
    if replay.source_copies() != plan.copies {
        return Err(invalid(
            "candidate copy placements disagree with the actual graph replay".to_string(),
        ));
    }
    let len = target.sequence.len() as isize;
    let mut hits: BTreeMap<String, Vec<PafHit>> = BTreeMap::new();
    for copy in &plan.copies {
        let span = usize::try_from(copy.end - copy.start)
            .map_err(|_| invalid("negative projected node span".to_string()))?;
        let start = if contract.alignment_reverse_complemented {
            (-(frame_offset as isize + copy.end)).rem_euclid(len)
        } else {
            (frame_offset as isize + copy.start).rem_euclid(len)
        } as usize;
        let reverse = copy.reverse ^ contract.alignment_reverse_complemented;
        hits.entry(copy.node.clone()).or_default().push(PafHit {
            query_name: copy.node.clone(),
            query_length: span,
            query_start: 0,
            query_end: span,
            strand: if reverse { '-' } else { '+' },
            target_start: start,
            target_end: start + span,
            matches: span,
            block_length: span,
            mapq: 0,
            substitution_matches: span,
            substitution_mismatches: 0,
            primary: true,
            exact_flanks: false,
        });
    }
    let mapping = Mapping {
        by_node: hits
            .into_iter()
            .map(|(node, placements)| {
                (
                    node,
                    NodeMapping {
                        selected: placements.first().cloned(),
                        accepted_hits: placements.clone(),
                        all_hits: placements,
                    },
                )
            })
            .collect(),
        reference_len: target.sequence.len(),
        command: format!(
            "{}; exact_graph_replay; evidence={}",
            evidence.method,
            paths.debug_dir.display()
        ),
        stderr: String::new(),
    };
    let classes = classify_nodes(graph, &mapping);
    let rows = graph.order.iter().map(|node| VerifiedNodeSourceRow {
        node: node.clone(), node_class: classes[node].clone(),
        source: if plan.alignment_cigar == "." { "authoritative_fasta_exact_graph_replay" } else { "authoritative_fasta_recorded_resolve_path" }.to_string(),
        length: plan.sequences[node].len(),
        parts: mapping.by_node[node].accepted_hits.iter().map(|hit| format!(
            "{}-{}:{}", hit.target_start + 1, hit.target_end, hit.strand,
        )).collect::<Vec<_>>().join(";"),
        notes: if plan.alignment_cigar == "." {
            "one-based unfolded intervals from exact graph replay; no whole-path alignment CIGAR was inferred; all shared copies agree"
        } else {
            "zero-based projection converted to one-based unfolded intervals; placements derived from the whole-path alignment; all shared copies agree"
        }.to_string(),
    }).collect();
    let changed = plan
        .edge_overlaps
        .iter()
        .filter(|edge| edge.old_overlap != edge.new_overlap)
        .collect::<Vec<_>>();
    let mut summary = vec![
        ("topology_contract", "PASS_checked_compaction".to_string()),
        (
            "recorded_path_replay",
            "PASS_exact_authoritative_sequence".to_string(),
        ),
        ("all_edge_sequence_overlaps", "PASS".to_string()),
        (
            "shared_node_sequences",
            "PASS_all_recorded_occurrences".to_string(),
        ),
        ("circle_evidence", contract.circular_evidence.clone()),
        (
            "restored_checked_links",
            contract.restored_links.len().to_string(),
        ),
        (
            "checked_compaction_nodes",
            contract.lineage.len().to_string(),
        ),
        (
            "recorded_path_occurrences",
            plan.copies.len().to_string(),
        ),
        (
            "source_path_length",
            evidence.source_path_length.map_or_else(|| ".".to_string(), |length| length.to_string()),
        ),
        (
            "source_reference_placement_fallback",
            contract.reference_placement_fallback.to_string(),
        ),
        (
            "source_alignment_reverse_complemented",
            contract.alignment_reverse_complemented.to_string(),
        ),
        (
            "source_alignment_rotation_step",
            contract.alignment_rotation_step.to_string(),
        ),
        ("target_frame_offset", frame_offset.to_string()),
        ("path_alignment_method", alignment_method),
        ("accepted_path_alignments", evidence.accepted_alignments.to_string()),
        ("multiple_valid_assignments", evidence.multiple_valid_assignments.to_string()),
        ("candidate_generation_method", evidence.method.clone()),
        ("path_alignment_cigar", plan.alignment_cigar.clone()),
        ("alignment_candidates_examined", plan.alignment_candidates_examined.to_string()),
        (
            "path_assignment_selection",
            if plan.alignment_cigar == "." {
                "supplied_node_candidate_passed_exact_graph_replay"
            } else {
                "minimum_changed_overlaps,overlap_delta,changed_nodes,node_length_delta,changed_sequence_span;exact_content_tiebreak"
            }.to_string(),
        ),
        (
            "selected_assignment_preference",
            format!("{},{},{},{},{}", selected_preference.changed_overlaps,
                selected_preference.overlap_delta, selected_preference.changed_nodes,
                selected_preference.node_length_delta, selected_preference.changed_sequence_span),
        ),
        (
            "boundary_assignments_explored",
            plan.explored_assignments.to_string(),
        ),
        ("changed_overlap_count", changed.len().to_string()),
        (
            "changed_overlaps",
            changed
                .iter()
                .map(|edge| {
                    format!(
                        "{}:{}M>{}M",
                        edge.edge_index, edge.old_overlap, edge.new_overlap
                    )
                })
                .collect::<Vec<_>>()
                .join(";"),
        ),
        (
            "linear_tiling_role",
            "diagnostic_only; graph overlaps may cover a base more than once".to_string(),
        ),
    ];
    // Reconstruction provenance (appended AFTER the original 25 recorded-path fields so
    // the original schema/order is untouched). generator1 = whole-path recorded-walk
    // alignment; generator2 = independent per-node projection with exact replay. For
    // generator1 we additionally record the scenario: ordinary, A (shared repeat-node
    // boundary repaired by homologous-core trimming) and/or B (tolerated source overlap
    // end-sequence conflict). A true source-graph structural error (C) cannot reach a
    // success report; it is flagged on failure instead.
    let is_generator1 = evidence.method == "recorded_resolve_path_alignment";
    let homologous_core_applied = !plan.shared_core_trims.is_empty();
    let source_conflict_present = !evidence.source_overlap_conflicts.is_empty();
    let mut scenario_parts = Vec::new();
    if homologous_core_applied {
        scenario_parts.push("A_shared_repeat_boundary");
    }
    if source_conflict_present {
        scenario_parts.push("B_source_overlap_conflict");
    }
    let reconstruction_scenario = if !is_generator1 {
        ".".to_string()
    } else if scenario_parts.is_empty() {
        "ordinary".to_string()
    } else {
        scenario_parts.join("+")
    };
    let homologous_core_trims = if homologous_core_applied {
        plan.shared_core_trims
            .iter()
            .map(|trim| {
                // source_length is the ORIGINAL repeat-node span (its minimap2/blastn
                // self-alignment boundary); the retained 0-based interval/new length is
                // the narrower homologous core. Keep both for future repeat up/downgrade.
                format!(
                    "{}:source_length={}:retained_0based={}-{}:retained_length={}",
                    trim.node, trim.source_length, trim.source_start, trim.source_end, trim.new_length
                )
            })
            .collect::<Vec<_>>()
            .join(";")
    } else {
        ".".to_string()
    };
    summary.push((
        "reconstruction_generator",
        if is_generator1 {
            "generator1_recorded_path"
        } else {
            "generator2_node_projection"
        }
        .to_string(),
    ));
    summary.push(("reconstruction_scenario", reconstruction_scenario));
    summary.push((
        "homologous_core_applied",
        homologous_core_applied.to_string(),
    ));
    summary.push(("homologous_core_trims", homologous_core_trims));
    summary.push((
        "homologous_core_note",
        if homologous_core_applied {
            "retained homologous core may be narrower than the minimap2/blastn self-alignment repeat boundary; source_length retained for future repeat upgrade/downgrade reconciliation"
        } else {
            "."
        }
        .to_string(),
    ));
    graph.links = updated_links;
    // Internal-only diagnostics (homologous shared-core trims and tolerated source
    // overlap end-sequence conflicts) are deliberately kept OUT of the run report so
    // the final report schema stays identical to the original rebuild; overlap changes
    // still surface through the original `changed_overlaps` fields. They are audited
    // here, in the debug directory, for traceability.
    let mut internal_audit = String::from("kind\tdetail\n");
    for trim in &plan.shared_core_trims {
        internal_audit.push_str(&format!(
            "shared_core_trim\t{}:source_length={}:source_interval_0based={}-{}:new_length={}\n",
            trim.node, trim.source_length, trim.source_start, trim.source_end, trim.new_length
        ));
    }
    for conflict in &evidence.source_overlap_conflicts {
        internal_audit.push_str(&format!(
            "source_overlap_conflict\tedge{}:{}>{}:declared_{}M:mismatch_{}\n",
            conflict.edge_index, conflict.from, conflict.to,
            conflict.declared_overlap, conflict.mismatch_count
        ));
    }
    fs::create_dir_all(&paths.debug_dir)?;
    fs::write(
        paths.debug_dir.join("recorded_path_internal_audit.tsv"),
        internal_audit,
    )?;
    Ok((
        SequencePlan {
            sequences: plan.sequences,
            mapping,
            node_source_rows: rows,
            method: evidence.method,
        },
        classes,
        summary,
    ))
}

#[cfg(test)]
mod tests_recorded_path {
    use super::*;

    fn candidate_fixture(
        nodes: &[(&str, &str)],
        edges: &[(&str, char, &str, char, usize)],
        walk: &[(&str, char)],
        circular: bool,
    ) -> (Gfa, ResolveContract) {
        let source_sequences = nodes
            .iter()
            .map(|(name, sequence)| (name.to_string(), sequence.to_string()))
            .collect::<BTreeMap<_, _>>();
        let topology_links = edges
            .iter()
            .map(|(a, ao, b, bo, n)| EvidenceLink {
                from: a.to_string(),
                from_orientation: *ao,
                to: b.to_string(),
                to_orientation: *bo,
                overlap: format!("{n}M"),
                tags: Vec::new(),
            })
            .collect::<Vec<_>>();
        let graph = Gfa {
            headers: Vec::new(),
            segments: source_sequences
                .iter()
                .map(|(name, sequence)| {
                    (
                        name.clone(),
                        Segment {
                            name: name.clone(),
                            sequence: sequence.clone(),
                            tags: Vec::new(),
                        },
                    )
                })
                .collect(),
            order: nodes.iter().map(|(name, _)| name.to_string()).collect(),
            links: topology_links
                .iter()
                .map(|edge| Link {
                    from_name: edge.from.clone(),
                    from_orient: edge.from_orientation,
                    to_name: edge.to.clone(),
                    to_orient: edge.to_orientation,
                    overlap: edge.overlap.clone(),
                    tags: Vec::new(),
                })
                .collect(),
            other_lines: Vec::new(),
        };
        let contract = ResolveContract {
            source_sequences,
            resolved_links: topology_links.clone(),
            topology_links,
            restored_links: Vec::new(),
            lineage: BTreeMap::new(),
            walk: walk
                .iter()
                .enumerate()
                .map(|(i, (node, sign))| WalkStep {
                    source_node: node.to_string(),
                    instance_node: format!("{node}_occurrence_{i}"),
                    orientation: *sign,
                })
                .collect(),
            circular,
            circular_evidence: "test_explicit_graph_evidence".to_string(),
            reference_placement_fallback: false,
            alignment_reverse_complemented: false,
            alignment_rotation_step: 0,
            path_issues: Vec::new(),
        };
        (graph, contract)
    }

    fn accept_candidate(
        graph: &mut Gfa,
        contract: &ResolveContract,
        sequences: BTreeMap<String, String>,
        target: &str,
    ) -> Result<BuildResult, OrgraftError> {
        accept_node_sequences(
            graph,
            &FastaRecord {
                header: "test".to_string(),
                id: "test".to_string(),
                sequence: target.to_string(),
            },
            contract,
            sequences,
            // Tests must never materialise output inside the repository working
            // tree; route the candidate-verification scratch dir to the OS temp
            // location with a per-process, per-call unique name.
            &OutputPaths::new(
                &std::env::temp_dir().join(format!(
                    "orgraft_candidate_verification_{}_{}",
                    std::process::id(),
                    std::time::SystemTime::now()
                        .duration_since(std::time::UNIX_EPOCH)
                        .unwrap_or_default()
                        .as_nanos()
                )),
                "subgraph_001",
            ),
            "node_projection_with_exact_graph_replay",
        )
    }

    struct JointFixture {
        root: PathBuf,
        graph: Gfa,
        contract: ResolveContract,
        sequences: BTreeMap<String, String>,
        target: String,
    }

    impl Drop for JointFixture {
        fn drop(&mut self) { let _ = fs::remove_dir_all(&self.root); }
    }

    impl JointFixture {
        fn new(reverse_arms: bool) -> Self {
            let root = std::env::temp_dir().join(format!("orgraft_joint_ir_{}_{}",
                std::process::id(), SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos()));
            fs::create_dir(&root).unwrap();
            let a = "A".repeat(2400);
            let b = "C".repeat(2400);
            let (graph, contract) = candidate_fixture(
                &[("a", &a), ("r", "ACGTAC"), ("b", &b)],
                &[("a", '+', "r", '+', 0), ("r", '+', "b", '+', 0),
                  ("b", '+', "r", '-', 0), ("r", '-', "a", '+', 0)],
                &[("a", '+'), ("r", '+'), ("b", '+'), ("r", '-')], true);
            let repeat = "ACGTACGA".to_string();
            let target = format!("{a}{repeat}{b}{}", reverse_complement(&repeat));
            let sequences = [("a".into(), a), ("r".into(), repeat), ("b".into(), b)].into();
            fs::write(root.join("checked.gfa"), "original checked input\n").unwrap();
            fs::write(root.join("resolve_details.tsv"), "original recorded resolve input\n").unwrap();
            fs::write(root.join("candidate.fasta"), format!(">subgraph_001\n{target}\n")).unwrap();
            let fixture = Self { root, graph, contract, sequences, target };
            let mut trial = fixture.sequences.clone();
            if reverse_arms { for name in ["a", "b"] {
                trial.insert(name.into(), reverse_complement(&trial[name]));
            }}
            let mut text = String::new();
            for (name, sequence) in trial { writeln!(text, "S\t{name}\t{sequence}").unwrap(); }
            for link in &fixture.graph.links {
                writeln!(text, "L\t{}\t{}\t{}\t{}\t{}", link.from_name, link.from_orient,
                    link.to_name, link.to_orient, link.overlap).unwrap();
            }
            fs::write(fixture.root.join("candidate.gfa"), text).unwrap();
            fixture.rebind();
            fixture
        }

        fn rebind(&self) {
            let mut text = "key\tvalue\n".to_string();
            for (key, file) in [("checked_gfa_sha256", "checked.gfa"),
                ("resolve_details_sha256", "resolve_details.tsv"),
                ("candidate_fasta_sha256", "candidate.fasta"),
                ("candidate_gfa_sha256", "candidate.gfa")] {
                writeln!(text, "{key}\t{}", file_sha(&self.root.join(file)).unwrap()).unwrap();
            }
            text.push_str("status\tselected_for_ordinary_validation\nR2_path_links\tverified_zero_overlap\nopposite_pair_count\t1\nfixed_read_pool\t10\nmain_ALT\t6\nother_complete_ALT\t0\nremainder\t4\nstrict_formula\t6 > 4\nrepeat_node\tr\n");
            fs::write(self.root.join("joint_ir.tsv"), text).unwrap();
        }

        fn validate(&self) -> Result<(BTreeMap<String, String>, bool), OrgraftError> {
            joint_ir_trial_sequences(&self.root.join("candidate.gfa"), &self.root.join("joint_ir.tsv"),
                &self.root.join("checked.gfa"), &self.root.join("resolve_details.tsv"),
                &self.root.join("candidate.fasta"), &self.contract)
        }
    }

    #[test]
    fn joint_ir_keeps_shared_repeat_and_replays_the_original_checked_path() {
        for reverse in [false, true] {
            let fixture = JointFixture::new(reverse);
            let (sequences, normalized) = fixture.validate().unwrap();
            assert_eq!(sequences, fixture.sequences);
            assert_eq!(normalized, reverse);
            let mut graph = fixture.graph.clone();
            let (_, _, summary) = accept_node_sequences(&mut graph,
                &FastaRecord { header: "test".into(), id: "test".into(), sequence: fixture.target.clone() },
                &fixture.contract, sequences.clone(), &OutputPaths::new(&fixture.root.join("rebuild"), "subgraph_001"),
                "joint_ir_hash_bound_original_checked_path").unwrap();
            assert_eq!(graph.segments.len(), 3);
            assert!(summary.iter().any(|(key,value)| *key == "candidate_generation_method"
                && value == "joint_ir_hash_bound_original_checked_path"));
            let mut wrong = sequences;
            wrong.get_mut("r").unwrap().push('A');
            assert!(accept_candidate(&mut fixture.graph.clone(), &fixture.contract, wrong, &fixture.target).is_err());
        }
    }

    #[test]
    fn joint_ir_rejects_changed_hashes_and_duplicate_evidence_keys() {
        for file in ["checked.gfa", "resolve_details.tsv", "candidate.fasta", "candidate.gfa"] {
            let fixture = JointFixture::new(false);
            OpenOptions::new().append(true).open(fixture.root.join(file)).unwrap().write_all(b"\n").unwrap();
            assert!(fixture.validate().unwrap_err().to_string().contains("input hash differs"));
        }
        let fixture = JointFixture::new(false);
        OpenOptions::new().append(true).open(fixture.root.join("joint_ir.tsv")).unwrap()
            .write_all(b"main_ALT\t6\n").unwrap();
        assert!(fixture.validate().unwrap_err().to_string().contains("duplicate event key"));
    }

    #[test]
    fn joint_ir_rejects_failed_fixed_pool_and_unselected_events() {
        for (old, new, expected) in [
            ("main_ALT\t6", "main_ALT\t4", "strict gate"),
            ("fixed_read_pool\t10", "fixed_read_pool\t11", "strict gate"),
            ("strict_formula\t6 > 4", "strict_formula\t6 >= 4", "strict gate"),
            ("status\tselected_for_ordinary_validation", "status\treview_required", "not selected"),
            ("opposite_pair_count\t1", "opposite_pair_count\t2", "not selected"),
        ] {
            let fixture = JointFixture::new(false);
            let path = fixture.root.join("joint_ir.tsv");
            let text = fs::read_to_string(&path).unwrap();
            fs::write(path, text.replace(old,new)).unwrap();
            assert!(fixture.validate().unwrap_err().to_string().contains(expected));
        }
    }

    #[test]
    fn joint_ir_rejects_rebound_topology_changes_and_disagreeing_arm_orientation() {
        let fixture = JointFixture::new(false);
        let path = fixture.root.join("candidate.gfa");
        let text = fs::read_to_string(&path).unwrap();
        fs::write(&path, text.replace("L\ta\t+\tr\t+\t0M", "L\ta\t+\tr\t+\t1M")).unwrap();
        fixture.rebind();
        assert!(fixture.validate().unwrap_err().to_string().contains("connections or overlaps"));
        fs::write(&path, text.replace(&format!("S\ta\t{}", "A".repeat(2400)),
            &format!("S\ta\t{}", "T".repeat(2400)))).unwrap();
        fixture.rebind();
        assert!(fixture.validate().unwrap_err().to_string().contains("disagree on orientation"));
    }

    #[test]
    fn joint_ir_cli_requires_complete_original_source_and_fresh_output() {
        let base = ["--edited-gfa", "original.gfa", "--polished-fasta", "candidate.fasta",
            "--checked-gfa", "checked.gfa", "--resolve-details", "resolve.tsv",
            "--joint-ir-trial-gfa", "candidate.gfa", "--joint-ir-evidence", "joint_ir.tsv"]
            .map(str::to_string).to_vec();
        assert!(RebuildOptions::from_args(&base).is_ok());
        for flag in ["--force", "--overlap-only"] {
            let mut args = base.clone(); args.push(flag.into());
            assert!(RebuildOptions::from_args(&args).is_err());
        }
        for key in ["--joint-ir-evidence", "--joint-ir-trial-gfa", "--checked-gfa", "--resolve-details"] {
            let mut args = base.clone(); let index = args.iter().position(|v| v == key).unwrap();
            args.drain(index..index+2);
            assert!(RebuildOptions::from_args(&args).is_err());
        }
    }

    #[test]
    fn node_candidate_requires_exact_replay_and_does_not_invent_a_cigar() {
        let (mut graph, contract) = candidate_fixture(
            &[("a", "AAC"), ("b", "CGG")],
            &[("a", '+', "b", '+', 1)],
            &[("a", '+'), ("b", '+')],
            false,
        );
        let candidate = BTreeMap::from([
            ("a".to_string(), "AATC".to_string()),
            ("b".to_string(), "CGG".to_string()),
        ]);
        let (plan, _, rows) = accept_candidate(&mut graph, &contract, candidate, "AATCGG").unwrap();
        let summary = rows.into_iter().collect::<BTreeMap<_, _>>();
        assert_eq!(plan.method, "node_projection_with_exact_graph_replay");
        assert_eq!(summary["path_alignment_cigar"], ".");
        assert_eq!(summary["alignment_candidates_examined"], "0");
        assert_eq!(summary["accepted_path_alignments"], "0");
        assert_eq!(plan.mapping.by_node["a"].accepted_hits[0].target_end, 4);
        assert_eq!(plan.mapping.by_node["b"].accepted_hits[0].target_start, 3);
    }

    #[test]
    fn candidate_checks_an_unused_edge_even_when_the_selected_walk_matches() {
        let (mut graph, contract) = candidate_fixture(
            &[("a", "AAC"), ("b", "CGG")],
            &[("a", '+', "b", '+', 1), ("b", '+', "a", '+', 1)],
            &[("a", '+'), ("b", '+')],
            false,
        );
        let error = accept_candidate(
            &mut graph,
            &contract,
            contract.source_sequences.clone(),
            "AACGG",
        )
        .unwrap_err();
        assert!(error
            .to_string()
            .contains("candidate edge overlap conflicts"));
    }

    #[test]
    fn circular_candidate_removes_one_closing_overlap_and_maps_the_full_node() {
        let (mut graph, contract) = candidate_fixture(
            &[("a", "ACGTAC")],
            &[("a", '+', "a", '+', 2)],
            &[("a", '+')],
            true,
        );
        let (plan, _, rows) = accept_candidate(
            &mut graph,
            &contract,
            contract.source_sequences.clone(),
            "ACGT",
        )
        .unwrap();
        let summary = rows.into_iter().collect::<BTreeMap<_, _>>();
        assert_eq!(summary["target_frame_offset"], "2");
        let hit = &plan.mapping.by_node["a"].accepted_hits[0];
        assert_eq!(
            (hit.target_start, hit.target_end, hit.query_length),
            (0, 6, 6)
        );
        assert!(accept_candidate(
            &mut graph,
            &contract,
            contract.source_sequences.clone(),
            "ACGTAC"
        )
        .is_err());
    }

    #[test]
    fn linear_candidate_does_not_gain_rotation_or_unrecorded_reverse_complement() {
        let (mut graph, mut contract) =
            candidate_fixture(&[("a", "AACG")], &[], &[("a", '+')], false);
        assert!(accept_candidate(
            &mut graph,
            &contract,
            contract.source_sequences.clone(),
            "ACGA"
        )
        .is_err());
        assert!(accept_candidate(
            &mut graph,
            &contract,
            contract.source_sequences.clone(),
            "CGTT"
        )
        .is_err());
        contract.alignment_reverse_complemented = true;
        let (plan, _, _) = accept_candidate(
            &mut graph,
            &contract,
            contract.source_sequences.clone(),
            "CGTT",
        )
        .unwrap();
        assert_eq!(plan.mapping.by_node["a"].accepted_hits[0].strand, '-');
    }

    #[test]
    fn node_candidate_cannot_change_connections_or_split_a_shared_node() {
        let (mut graph, contract) = candidate_fixture(
            &[("r", "AC"), ("b", "GG")],
            &[("r", '+', "b", '+', 0), ("b", '+', "r", '+', 0)],
            &[("r", '+'), ("b", '+'), ("r", '+')],
            false,
        );
        let mut candidate = contract.source_sequences.clone();
        candidate.insert("r_copy".to_string(), "AT".to_string());
        assert!(accept_candidate(&mut graph, &contract, candidate, "ACGGAT").is_err());
        assert!(accept_candidate(
            &mut graph,
            &contract,
            contract.source_sequences.clone(),
            "ACGGAT"
        )
        .is_err());
        let (plan, _, _) = accept_candidate(
            &mut graph,
            &contract,
            contract.source_sequences.clone(),
            "ACGGAC",
        )
        .unwrap();
        assert_eq!(plan.mapping.by_node["r"].accepted_hits.len(), 2);
        graph.links[0].to_orient = '-';
        let error = accept_candidate(
            &mut graph,
            &contract,
            contract.source_sequences.clone(),
            "ACGGAC",
        )
        .unwrap_err();
        assert!(error.to_string().contains("directed connections differ"));
    }

    fn circle_assignment(
        source: &str,
        target: &str,
        cigar: &str,
        offset: usize,
    ) -> AcceptedAssignment {
        let nodes = BTreeMap::from([("node".to_string(), source.to_string())]);
        let edges = vec![PathEdge {
            from: "node".to_string(),
            from_reverse: false,
            to: "node".to_string(),
            to_reverse: false,
            overlap: 0,
        }];
        let walk = vec![PathStep {
            node: "node".to_string(),
            reverse: false,
        }];
        let prepared = prepare_path(&nodes, &edges, &walk, true).unwrap();
        let plan =
            project_path_sequence(&prepared, target, cigar, 10_000).unwrap();
        let replay =
            prepare_path(&plan.sequences, &edges, &walk, true).unwrap();
        assert_eq!(replay.source_spell, target);
        (
            plan,
            offset,
            target.to_string(),
            "test_complete_alignment".to_string(),
        )
    }

    #[test]
    fn two_valid_circle_origins_do_not_require_one_node_representation() {
        let source = BTreeMap::from([("node".to_string(), "AGTC".to_string())]);
        let after = circle_assignment("AGTC", "AGTCA", "4=1D", 7);
        let before = circle_assignment("AGTC", "AAGTC", "1D4=", 6);
        assert_ne!(after.0.sequences, before.0.sequences);
        assert_eq!(after.0.edge_overlaps, before.0.edge_overlaps);
        assert!((after.2.clone() + &after.2).contains(&before.2));
        assert_eq!(assignment_preference(&after.0, &source).changed_overlaps, 0);
        assert_eq!(
            assignment_preference(&before.0, &source).changed_overlaps,
            0
        );
        let winner = |first: &AcceptedAssignment, second: &AcceptedAssignment| {
            if compare_assignments(second, first, &source).is_lt() {
                second.1
            } else {
                first.1
            }
        };
        assert_eq!(winner(&after, &before), winner(&before, &after));
    }

    #[test]
    fn preference_retains_original_node_origin_when_only_circle_rotation_differs() {
        let source = BTreeMap::from([("node".to_string(), "AGTC".to_string())]);
        let nearby = circle_assignment("AGTC", "AGTCA", "4=1D", 7);
        // Both targets are the same circle and both plans have a complete exact replay.
        let shifted = circle_assignment("AGTC", "GTCAA", "1I3=2D", 8);
        assert!((nearby.2.clone() + &nearby.2).contains(&shifted.2));
        assert!(compare_assignments(&nearby, &shifted, &source).is_lt());
    }

    #[test]
    fn preference_counts_changed_overlap_before_sequence_span() {
        let a = AssignmentPreference {
            changed_overlaps: 0,
            overlap_delta: 0,
            changed_nodes: 2,
            node_length_delta: 3,
            changed_sequence_span: 20,
        };
        let b = AssignmentPreference {
            changed_overlaps: 1,
            overlap_delta: 1,
            changed_nodes: 1,
            node_length_delta: 1,
            changed_sequence_span: 1,
        };
        assert!(a < b);
    }

    #[test]
    fn common_terminal_sequence_is_not_counted_twice() {
        let source = BTreeMap::from([("node".to_string(), "AAAA".to_string())]);
        let longer = circle_assignment("AAAA", "AAAAA", "4=1D", 0);
        let shorter = circle_assignment("AAAA", "AAA", "3=1I", 0);
        assert_eq!(
            assignment_preference(&longer.0, &source).changed_sequence_span,
            1
        );
        assert_eq!(
            assignment_preference(&shorter.0, &source).changed_sequence_span,
            1
        );
    }
}



// =====================================================================
// Inlined module: rebuild_structure.rs (manual structure corrections; ENTRY DISABLED)
// =====================================================================
// Apply explicitly supported graph corrections before sequence reconstruction.
// Evidence provenance is checked here; biological support is supplied by the
// reviewed correction record, not inferred from FASTA similarity alone.


#[derive(Debug)]
pub(super) struct AppliedCorrections {
    pub record: String,
    pub count: usize,
    pub removed_links: usize,
    pub added_links: usize,
    pub reversed_nodes: usize,
    pub removed_nodes: usize,
    /// Whole-node sequence replacements (node_sequence); internal bookkeeping only,
    /// not emitted as a new report key so the final report schema stays unchanged.
    #[allow(dead_code)]
    pub replaced_nodes: usize,
    /// Linear non-repeat merged nodes split into child nodes (split_node).
    #[allow(dead_code)]
    pub split_nodes: usize,
    pub old_walk: String,
    pub new_walk: String,
}

fn require(value: bool, message: impl Into<String>) -> Result<(), OrgraftError> {
    if value {
        Ok(())
    } else {
        Err(invalid(message))
    }
}
fn file_sha(path: &Path) -> Result<String, OrgraftError> {
    Ok(sha256_hex(&fs::read(path)?))
}
/// Read a small hash-bound FASTA of corrected node sequences into an uppercase map.
fn read_corrected_fasta(path: &Path) -> Result<BTreeMap<String, String>, OrgraftError> {
    let text = fs::read_to_string(path)?;
    let mut out: BTreeMap<String, String> = BTreeMap::new();
    let mut current: Option<String> = None;
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            let name = header
                .split_whitespace()
                .next()
                .ok_or_else(|| invalid("corrected FASTA has an empty header"))?
                .to_string();
            require(out.insert(name.clone(), String::new()).is_none(),
                "corrected FASTA has duplicate record")?;
            current = Some(name);
        } else if let Some(name) = &current {
            out.get_mut(name).unwrap().push_str(line.trim());
        } else {
            return Err(invalid("corrected FASTA has sequence before any header"));
        }
    }
    for (name, seq) in out.iter_mut() {
        *seq = seq.to_ascii_uppercase();
        require(!seq.is_empty(), format!("corrected FASTA record {name} is empty"))?;
        require(
            seq.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')),
            format!("corrected FASTA record {name} has non-IUPAC-ACGTN bases"),
        )?;
    }
    Ok(out)
}
/// Lazily load (and re-verify the bound hash of) a hash-bound corrected-sequence FASTA.
fn corrected_records<'a>(
    cache: &'a mut BTreeMap<String, BTreeMap<String, String>>,
    evidence_id: &str,
    evidence_paths: &BTreeMap<String, std::path::PathBuf>,
    actual_evidence: &BTreeMap<String, String>,
) -> Result<&'a BTreeMap<String, String>, OrgraftError> {
    if !cache.contains_key(evidence_id) {
        let path = evidence_paths
            .get(evidence_id)
            .ok_or_else(|| invalid("node correction references undeclared evidence"))?;
        require(
            file_sha(path)? == actual_evidence[evidence_id],
            "corrected FASTA hash changed at load time",
        )?;
        cache.insert(evidence_id.to_string(), read_corrected_fasta(path)?);
    }
    Ok(cache.get(evidence_id).unwrap())
}
fn endpoint(value: &str) -> Result<(String, char), OrgraftError> {
    require(value.len() > 1, "empty correction endpoint")?;
    let (node, sign) = value.split_at(value.len() - 1);
    require(
        matches!(sign, "+" | "-"),
        "invalid correction endpoint orientation",
    )?;
    Ok((node.to_string(), sign.chars().next().unwrap()))
}
fn key(link: &EvidenceLink) -> ((String, char), (String, char)) {
    let a = (link.from.clone(), link.from_orientation);
    let b = (link.to.clone(), flip(link.to_orientation));
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}
fn walk_text(walk: &[WalkStep]) -> String {
    walk.iter()
        .map(|w| format!("{}{}", w.instance_node, w.orientation))
        .collect::<Vec<_>>()
        .join(",")
}
fn branch_nodes(links: &[EvidenceLink]) -> BTreeSet<String> {
    let mut unique = BTreeSet::new();
    let mut degree = BTreeMap::<String, usize>::new();
    for edge in links {
        if unique.insert(key(edge)) {
            *degree.entry(edge.from.clone()).or_default() += 1;
            *degree.entry(edge.to.clone()).or_default() += 1;
        }
    }
    degree
        .into_iter()
        .filter_map(|(node, n)| if n >= 3 { Some(node) } else { None })
        .collect()
}

pub(super) fn apply(
    record_path: &Path,
    contract: &mut ResolveContract,
    bindings: &[(&str, &Path)],
) -> Result<AppliedCorrections, OrgraftError> {
    let text = fs::read_to_string(record_path)?;
    let hashes = bindings
        .iter()
        .map(|(role, path)| Ok((role.to_string(), file_sha(path)?)))
        .collect::<Result<BTreeMap<_, _>, OrgraftError>>()?;
    let mut evidence = BTreeMap::new();
    for line in text.lines().skip(1) {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.first() == Some(&"evidence") {
            require(fields.len() == 6, "invalid structure evidence row")?;
            require(
                Path::new(fields[3]).is_absolute(),
                "structure evidence requires an absolute source path",
            )?;
            evidence.insert(fields[1].to_string(), file_sha(Path::new(fields[3]))?);
        }
    }
    apply_text(&text, contract, &hashes, &evidence)
}

fn apply_text(
    text: &str,
    contract: &mut ResolveContract,
    actual_inputs: &BTreeMap<String, String>,
    actual_evidence: &BTreeMap<String, String>,
) -> Result<AppliedCorrections, OrgraftError> {
    let mut lines = text.lines();
    require(
        lines.next() == Some("kind\ttarget\tbefore\tafter\tevidence\tnote"),
        "invalid structure correction header",
    )?;
    let rows = lines
        .filter(|s| !s.trim().is_empty())
        .map(|line| line.split('\t').collect::<Vec<_>>())
        .collect::<Vec<_>>();
    require(
        rows.iter().all(|r| r.len() == 6),
        "structure correction rows require six columns",
    )?;
    let mut inputs = BTreeMap::new();
    let mut evidence = BTreeSet::new();
    let mut evidence_paths = BTreeMap::new();
    let mut actions = Vec::new();
    for row in &rows {
        require(
            !row[5].trim().is_empty() && row[5] != ".",
            "structure correction rows require a reason or source description",
        )?;
        match row[0] {
            "input" => {
                require(
                    inputs
                        .insert(row[1].to_string(), row[2].to_string())
                        .is_none(),
                    "duplicate correction input",
                )?;
                require(
                    row[3] == "." && row[4] == ".",
                    "invalid correction input row",
                )?;
            }
            "evidence" => {
                require(
                    evidence.insert(row[1].to_string()),
                    "duplicate structure evidence identifier",
                )?;
                require(
                    actual_evidence.get(row[1]).map(String::as_str) == Some(row[2]),
                    format!("structure evidence hash changed: {}", row[1]),
                )?;
                require(
                    Path::new(row[3]).is_absolute(),
                    "structure evidence corrected FASTA requires an absolute path",
                )?;
                evidence_paths.insert(row[1].to_string(), Path::new(row[3]).to_path_buf());
                require(
                    row[4] == ".",
                    "evidence records cannot refer to another evidence record",
                )?;
            }
            "link" | "reverse_node" | "remove_node" | "walk" | "circular"
            | "node_sequence" | "split_node" => actions.push(row),
            _ => {
                return Err(invalid(format!(
                    "unsupported structure correction kind: {}",
                    row[0]
                )))
            }
        }
    }
    let required = BTreeSet::from([
        "checked_gfa",
        "edited_gfa",
        "resolve_details",
        "polished_fasta",
    ]);
    require(
        inputs.keys().map(String::as_str).collect::<BTreeSet<_>>() == required,
        "structure corrections must bind all four graph/path/final inputs",
    )?;
    require(
        &inputs == actual_inputs,
        "structure correction input hash changed",
    )?;
    require(
        !evidence.is_empty() && !actions.is_empty(),
        "structure correction requires explicit evidence and an action",
    )?;
    for row in &actions {
        require(
            row[4].split(',').all(|id| evidence.contains(id)),
            "structure action references missing evidence",
        )?;
    }
    // Work on a clone so a malformed or incomplete correction never partly applies.
    let mut proposed = contract.clone();
    let source_branches = branch_nodes(&proposed.topology_links);
    let old_walk = walk_text(&proposed.walk);
    let old_counts = proposed
        .walk
        .iter()
        .fold(BTreeMap::<String, usize>::new(), |mut m, w| {
            *m.entry(w.source_node.clone()).or_default() += 1;
            m
        });
    let no_aliases = proposed
        .walk
        .iter()
        .all(|w| w.source_node == w.instance_node);
    let mut changed = BTreeSet::new();
    let mut removed_nodes = BTreeSet::new();
    // Parent -> ordered child names produced by split_node, and every new child node.
    let mut split_children: BTreeMap<String, Vec<String>> = BTreeMap::new();
    let mut new_child_nodes: BTreeSet<String> = BTreeSet::new();
    let mut replaced_nodes = 0usize;
    let mut split_nodes = 0usize;
    // Lazily loaded, hash-bound corrected-sequence FASTA per evidence identifier.
    let mut corrected_cache: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
    let mut counts = [0usize; 4];
    // Nodes whose incident topology is changed by an explicit (hash-bound) link row or
    // by a split. Only such nodes are allowed to gain/lose branching membership, so a
    // repeat can never be silently reclassified without a reviewed structure action.
    let mut edge_touched: BTreeSet<String> = BTreeSet::new();
    // Resolve a walk/link token to its source base node. A trailing `_copyK` marks a
    // repeat instance only when the stripped prefix is itself a real graph node; a base
    // node may literally carry a `_copy1` name from resolve (e.g. 64_65_66_copy1).
    let source_base = |raw: &str, seqs: &BTreeMap<String, String>| -> String {
        if let Some((base, suffix)) = raw.rsplit_once("_copy") {
            if !suffix.is_empty()
                && suffix.bytes().all(|b| b.is_ascii_digit())
                && seqs.contains_key(base)
            {
                return base.to_string();
            }
        }
        raw.to_string()
    };
    for row in &actions {
        require(
            changed.insert((row[0], row[1])),
            "duplicate structure action on the same object",
        )?;
        match row[0] {
            "link" => {
                let ends = row[1].split('>').collect::<Vec<_>>();
                require(
                    ends.len() == 2,
                    "link correction requires from+>to+ endpoints",
                )?;
                let (from, fo) = endpoint(ends[0])?;
                let (to, to_o) = endpoint(ends[1])?;
                // Link rows operate on source-node topology. Editing an instance alias
                // (_copyK whose stripped base is a real node) is still forbidden; base
                // source nodes are allowed even when the recorded walk uses aliases, and
                // a base node may itself carry a `_copy1`-style name from resolve.
                require(
                    source_base(&from, &proposed.source_sequences) == from
                        && source_base(&to, &proposed.source_sequences) == to,
                    "link corrections on resolved repeat aliases are not supported; repeat reclassification is outside this mode",
                )?;
                edge_touched.insert(from.clone());
                edge_touched.insert(to.clone());
                let link = EvidenceLink {
                    from,
                    from_orientation: fo,
                    to,
                    to_orientation: to_o,
                    overlap: row[3].to_string(),
                    tags: Vec::new(),
                };
                let wanted = key(&link);
                let existing = proposed
                    .topology_links
                    .iter()
                    .filter(|e| key(e) == wanted)
                    .collect::<Vec<_>>();
                if row[2] == "." {
                    require(
                        existing.is_empty(),
                        "structure link addition already exists",
                    )?;
                } else {
                    require(
                        existing.len() == 1 && existing[0].overlap == row[2],
                        "structure link before value does not match source",
                    )?;
                }
                require(row[2] != row[3], "structure link action has no change")?;
                if row[3] != "." {
                    link.overlap_bases()?;
                }
                if row[2] != "." {
                    counts[0] += 1;
                }
                if row[3] != "." {
                    counts[1] += 1;
                }
                proposed.topology_links.retain(|e| key(e) != wanted);
                if row[3] != "." {
                    proposed.topology_links.push(link);
                }
            }
            "reverse_node" => {
                let source = proposed
                    .source_sequences
                    .get(row[1])
                    .ok_or_else(|| invalid("structure node is absent"))?;
                require(
                    sha256_hex(source.as_bytes()) == row[2],
                    "structure node before sequence hash changed",
                )?;
                let reversed = reverse_complement(source);
                require(
                    sha256_hex(reversed.as_bytes()) == row[3],
                    "structure node after sequence hash does not match reverse complement",
                )?;
                require(
                    source != &reversed,
                    "structure node reversal has no sequence change",
                )?;
                proposed
                    .source_sequences
                    .insert(row[1].to_string(), reversed);
                counts[2] += 1;
            }
            "remove_node" => {
                require(
                    no_aliases
                        && old_counts.get(row[1]) == Some(&1)
                        && !source_branches.contains(row[1]),
                    "only a non-repeat node with explicit evidence can be removed",
                )?;
                let source = proposed
                    .source_sequences
                    .get(row[1])
                    .ok_or_else(|| invalid("removed structure node is absent"))?;
                require(
                    sha256_hex(source.as_bytes()) == row[2] && row[3] == ".",
                    "removed node sequence hash changed",
                )?;
                proposed.source_sequences.remove(row[1]);
                removed_nodes.insert(row[1].to_string());
                counts[3] += 1;
            }
            "walk" => {
                require(row[1] == "path", "walk correction target must be path")?;
                // Repeat aliases are allowed here: non-repeat copy counts are still
                // preserved below, branching membership is held fixed, and every
                // corrected transition must have a supporting source topology edge when
                // the instance-level resolved links are rebuilt.
                require(
                    row[2] == old_walk,
                    "structure walk before value differs from recorded resolve walk",
                )?;
                proposed.walk = row[3]
                    .split(',')
                    .map(|part| {
                        let (raw, orientation) = endpoint(part)?;
                        // A repeated copy is written "<base>_copyK"; its source node is
                        // the base node only when that base actually exists, otherwise the
                        // token is a base node that carries a `_copyK`-style name itself.
                        let source = source_base(&raw, &proposed.source_sequences);
                        Ok(WalkStep {
                            source_node: source,
                            instance_node: raw,
                            orientation,
                        })
                    })
                    .collect::<Result<Vec<_>, OrgraftError>>()?;
            }
            "circular" => {
                require(row[1] == "path", "circle correction target must be path")?;
                require(
                    matches!(row[2], "true" | "false") && matches!(row[3], "true" | "false"),
                    "invalid structure circle value",
                )?;
                require(
                    row[2] == proposed.circular.to_string() && row[2] != row[3],
                    "structure circle before value changed or action is empty",
                )?;
                proposed.circular = row[3] == "true";
                proposed.circular_evidence = "explicit_hash_bound_structure_correction".to_string();
            }
            "node_sequence" => {
                let node = row[1];
                require(
                    proposed.source_sequences.contains_key(node),
                    "node_sequence target node is missing",
                )?;
                let before = sha256_hex(proposed.source_sequences[node].as_bytes());
                require(row[2] == before, "node_sequence before hash changed")?;
                require(
                    evidence.contains(row[4]),
                    "node_sequence evidence identifier is not declared",
                )?;
                let records = corrected_records(
                    &mut corrected_cache,
                    row[4],
                    &evidence_paths,
                    actual_evidence,
                )?;
                let seq = records
                    .get(node)
                    .ok_or_else(|| invalid("corrected FASTA lacks the node_sequence record"))?;
                let after = sha256_hex(seq.as_bytes());
                require(row[3] == after, "node_sequence after hash does not match FASTA")?;
                require(after != before, "node_sequence action does not change the node")?;
                proposed.source_sequences.insert(node.to_string(), seq.clone());
                replaced_nodes += 1;
            }
            "split_node" => {
                let parts = row[1].split('>').collect::<Vec<_>>();
                require(
                    parts.len() == 2,
                    "split_node target must be parent>child1,child2",
                )?;
                let parent = parts[0];
                let children = parts[1]
                    .split(',')
                    .map(str::to_string)
                    .collect::<Vec<_>>();
                require(children.len() >= 2, "split_node needs at least two children")?;
                require(
                    proposed.source_sequences.contains_key(parent),
                    "split_node parent is missing",
                )?;
                require(
                    !source_branches.contains(parent),
                    "cannot split a branching node; repeat reclassification is outside this mode",
                )?;
                require(
                    old_counts.get(parent) == Some(&1),
                    "cannot split a repeat (multi-copy) node",
                )?;
                let before = sha256_hex(proposed.source_sequences[parent].as_bytes());
                require(row[2] == before, "split_node before hash changed")?;
                let after_hashes = row[3].split(',').collect::<Vec<_>>();
                require(
                    after_hashes.len() == children.len(),
                    "split_node child/hash count mismatch",
                )?;
                require(
                    evidence.contains(row[4]),
                    "split_node evidence identifier is not declared",
                )?;
                let records = corrected_records(
                    &mut corrected_cache,
                    row[4],
                    &evidence_paths,
                    actual_evidence,
                )?;
                for (child, hash) in children.iter().zip(after_hashes) {
                    require(
                        !proposed.source_sequences.contains_key(child),
                        "split_node child already exists in the graph",
                    )?;
                    require(
                        new_child_nodes.insert(child.clone()),
                        "split_node child declared more than once",
                    )?;
                    let seq = records
                        .get(child)
                        .ok_or_else(|| invalid("corrected FASTA lacks a split child record"))?;
                    require(
                        hash == sha256_hex(seq.as_bytes()),
                        "split_node child after hash does not match FASTA",
                    )?;
                    proposed.source_sequences.insert(child.clone(), seq.clone());
                    if let Some(class) = proposed.lineage.get(parent) {
                        proposed.lineage.insert(child.clone(), class.clone());
                    }
                }
                // Parent sequence is removed; its incident links must be removed and the
                // children reconnected by accompanying link rows, and the walk row re-routed.
                proposed.source_sequences.remove(parent);
                removed_nodes.insert(parent.to_string());
                split_children.insert(parent.to_string(), children.clone());
                edge_touched.insert(parent.to_string());
                for child in &children {
                    edge_touched.insert(child.clone());
                }
                split_nodes += 1;
            }
            _ => unreachable!(),
        }
    }
    require(
        !proposed.source_sequences.is_empty(),
        "structure correction removes the whole graph",
    )?;
    for link in &proposed.topology_links {
        require(
            proposed.source_sequences.contains_key(&link.from)
                && proposed.source_sequences.contains_key(&link.to),
            "structure correction leaves an edge referencing a missing node",
        )?;
    }
    let mut expected_counts = old_counts;
    for node in removed_nodes {
        expected_counts.remove(&node);
    }
    // Split children are linear single-occurrence nodes in the corrected walk.
    for children in split_children.values() {
        for child in children {
            expected_counts.insert(child.clone(), 1);
        }
    }
    let new_counts = proposed
        .walk
        .iter()
        .fold(BTreeMap::<String, usize>::new(), |mut m, w| {
            *m.entry(w.source_node.clone()).or_default() += 1;
            m
        });
    // Non-branch nodes keep an exact traversal count and must not be dropped. An
    // existing branch (repeat) node may legitimately be re-traversed a different
    // number of times when the recorded walk under/over-visited a repeated copy —
    // its branching degree already represents those copies. Branching membership
    // itself is still required to be unchanged by the check further below.
    for (node, cnt) in &new_counts {
        require(
            proposed.source_sequences.contains_key(node),
            "structure correction walk references a missing node",
        )?;
        if !source_branches.contains(node) {
            require(
                expected_counts.get(node) == Some(cnt),
                "structure correction changes non-repeat node copy counts or omits a source node",
            )?;
        }
    }
    for node in expected_counts.keys() {
        if !source_branches.contains(node) {
            require(
                new_counts.contains_key(node),
                "structure correction omits a non-repeat source node from the walk",
            )?;
        }
    }
    // Branching membership is held fixed except at nodes whose topology is changed by an
    // explicit hash-bound link/split action. This forbids silent repeat reclassification
    // while permitting an evidence-bound split to dissolve the false branch created by a
    // mis-merge (e.g. a small repeat wrongly routing through the merged parent).
    let new_branches = branch_nodes(&proposed.topology_links);
    for node in source_branches.symmetric_difference(&new_branches) {
        require(
            edge_touched.contains(node),
            "structure correction changes branching-node membership without an explicit link/split action; repeat reclassification is outside this mode",
        )?;
    }
    proposed
        .lineage
        .retain(|node, _| proposed.source_sequences.contains_key(node));
    let retained_edges = proposed
        .topology_links
        .iter()
        .map(key)
        .collect::<BTreeSet<_>>();
    proposed
        .restored_links
        .retain(|edge| retained_edges.contains(&key(edge)));
    if no_aliases {
        proposed.resolved_links = proposed.topology_links.clone();
    } else {
        // With repeat aliases, resolved links are instance-level. Rebuild them from the
        // corrected walk, taking each transition's overlap from the corrected source
        // topology (which the link rows have made target-consistent).
        let mut source_overlap = BTreeMap::new();
        for edge in &proposed.topology_links {
            require(
                source_overlap.insert(key(edge), edge.overlap.clone()).is_none(),
                "duplicate source topology edge after structure correction",
            )?;
        }
        let mut resolved = Vec::new();
        let mut add_transition = |left: &WalkStep, right: &WalkStep|
         -> Result<(), OrgraftError> {
            let source_edge = EvidenceLink {
                from: left.source_node.clone(),
                from_orientation: left.orientation,
                to: right.source_node.clone(),
                to_orientation: right.orientation,
                overlap: String::new(),
                tags: Vec::new(),
            };
            let overlap = source_overlap.get(&key(&source_edge)).ok_or_else(|| {
                invalid(format!(
                    "corrected walk transition {}{} -> {}{} has no supporting topology edge",
                    left.instance_node, left.orientation, right.instance_node, right.orientation
                ))
            })?;
            resolved.push(EvidenceLink {
                from: left.instance_node.clone(),
                from_orientation: left.orientation,
                to: right.instance_node.clone(),
                to_orientation: right.orientation,
                overlap: overlap.clone(),
                tags: Vec::new(),
            });
            Ok(())
        };
        for pair in proposed.walk.windows(2) {
            add_transition(&pair[0], &pair[1])?;
        }
        if proposed.circular {
            add_transition(proposed.walk.last().unwrap(), &proposed.walk[0])?;
        }
        proposed.resolved_links = resolved;
    }
    proposed.path_issues.clear();
    for pair in proposed.walk.windows(2) {
        proposed.edge_between(&pair[0], &pair[1])?;
    }
    if proposed.circular {
        proposed.edge_between(proposed.walk.last().unwrap(), &proposed.walk[0])?;
    }
    let result = AppliedCorrections {
        record: text.to_string(),
        count: actions.len(),
        removed_links: counts[0],
        added_links: counts[1],
        reversed_nodes: counts[2],
        removed_nodes: counts[3],
        replaced_nodes,
        split_nodes,
        old_walk,
        new_walk: walk_text(&proposed.walk),
    };
    *contract = proposed;
    Ok(result)
}

#[cfg(test)]
mod tests_structure {
    use super::*;
    fn graph() -> ResolveContract {
        let sequences = BTreeMap::from([
            ("a".to_string(), "ACG".to_string()),
            ("b".to_string(), "TTA".to_string()),
            ("c".to_string(), "GGC".to_string()),
        ]);
        let edge = |a: &str, b: &str| EvidenceLink {
            from: a.to_string(),
            from_orientation: '+',
            to: b.to_string(),
            to_orientation: '+',
            overlap: "0M".to_string(),
            tags: Vec::new(),
        };
        let links = vec![edge("a", "b"), edge("b", "c"), edge("c", "a")];
        ResolveContract {
            source_sequences: sequences,
            topology_links: links.clone(),
            resolved_links: links,
            restored_links: Vec::new(),
            lineage: BTreeMap::new(),
            walk: ["a", "c", "b"]
                .iter()
                .map(|n| WalkStep {
                    source_node: n.to_string(),
                    instance_node: n.to_string(),
                    orientation: '+',
                })
                .collect(),
            circular: true,
            circular_evidence: "recorded".to_string(),
            reference_placement_fallback: true,
            alignment_reverse_complemented: false,
            alignment_rotation_step: 0,
            path_issues: vec!["source fallback omitted edges".to_string()],
        }
    }
    fn record(actions: &str) -> (String, BTreeMap<String, String>, BTreeMap<String, String>) {
        let inputs = [
            "checked_gfa",
            "edited_gfa",
            "resolve_details",
            "polished_fasta",
        ]
        .into_iter()
        .map(|s| (s.to_string(), sha256_hex(s.as_bytes())))
        .collect::<BTreeMap<_, _>>();
        let evidence =
            BTreeMap::from([("native_graph".to_string(), sha256_hex(b"native evidence"))]);
        let mut text = "kind\ttarget\tbefore\tafter\tevidence\tnote\n".to_string();
        for (name, hash) in &inputs {
            text.push_str(&format!("input\t{name}\t{hash}\t.\t.\tunchanged source\n"));
        }
        text.push_str(&format!(
            "evidence\tnative_graph\t{}\t/evidence/source\t.\tvalidated independent source\n",
            evidence["native_graph"]
        ));
        text.push_str(actions);
        (text, inputs, evidence)
    }
    fn rewire() -> String {
        let mut s = String::new();
        for (a, b) in [("a", "b"), ("b", "c"), ("c", "a")] {
            s.push_str(&format!(
                "link\t{a}+>{b}+\t0M\t.\tnative_graph\tremove unsupported junction\n"
            ));
        }
        for (a, b) in [("a", "c"), ("c", "b"), ("b", "a")] {
            s.push_str(&format!(
                "link\t{a}+>{b}+\t.\t0M\tnative_graph\trestore supported junction\n"
            ));
        }
        s
    }
    #[test]
    fn supported_connections_make_the_recorded_walk_legal_without_changing_sequence() {
        let mut g = graph();
        let before = g.source_sequences.clone();
        let (text, inputs, evidence) = record(&rewire());
        let report = apply_text(&text, &mut g, &inputs, &evidence).unwrap();
        assert_eq!(g.source_sequences, before);
        assert_eq!(report.removed_links, 3);
        assert_eq!(report.added_links, 3);
        assert_eq!(report.old_walk, report.new_walk);
        assert!(g.path_issues.is_empty());
        for p in g.walk.windows(2) {
            g.edge_between(&p[0], &p[1]).unwrap();
        }
    }
    #[test]
    fn changed_input_hash_rejects_without_mutating_the_graph() {
        let mut g = graph();
        let before = g.topology_links.clone();
        let (text, mut inputs, evidence) = record(&rewire());
        inputs.insert("polished_fasta".to_string(), sha256_hex(b"another target"));
        assert!(apply_text(&text, &mut g, &inputs, &evidence).is_err());
        assert_eq!(g.topology_links, before);
    }
    #[test]
    fn unavailable_support_does_not_authorize_graph_changes() {
        let mut g = graph();
        let (text, inputs, _) = record(&rewire());
        assert!(apply_text(&text, &mut g, &inputs, &BTreeMap::new()).is_err());
    }
    #[test]
    fn incomplete_link_correction_is_rejected_atomically() {
        let mut g = graph();
        let before = g.topology_links.clone();
        let (text, inputs, evidence) =
            record("link\ta+>b+\t0M\t.\tnative_graph\tremove unsupported junction\n");
        assert!(apply_text(&text, &mut g, &inputs, &evidence).is_err());
        assert_eq!(g.topology_links, before);
    }
    #[test]
    fn explicit_node_reversal_changes_sequence_orientation_and_keeps_endpoints() {
        let mut g = graph();
        g.walk.swap(1, 2);
        g.path_issues.clear();
        let links = g.topology_links.clone();
        let old = sha256_hex(g.source_sequences["b"].as_bytes());
        let new = sha256_hex(b"TAA");
        let (text, inputs, evidence) = record(&format!(
            "reverse_node\tb\t{old}\t{new}\tnative_graph\tvalidated inversion of node sequence\n"
        ));
        let report = apply_text(&text, &mut g, &inputs, &evidence).unwrap();
        assert_eq!(g.source_sequences["b"], "TAA");
        assert_eq!(g.topology_links, links);
        assert_eq!(report.reversed_nodes, 1);
    }
    #[test]
    fn changing_path_copy_number_is_not_a_structure_correction() {
        let mut g = graph();
        let (text, inputs, evidence) =
            record("walk\tpath\ta+,c+,b+\ta+,b+,c+,b+\tnative_graph\tincorrect repeat addition\n");
        assert!(apply_text(&text, &mut g, &inputs, &evidence)
            .unwrap_err()
            .to_string()
            .contains("copy counts"));
    }
    #[test]
    fn removal_requires_incident_link_changes_and_one_source_occurrence() {
        let mut g = graph();
        let old = sha256_hex(g.source_sequences["b"].as_bytes());
        let (text, inputs, evidence) = record(&format!(
            "remove_node\tb\t{old}\t.\tnative_graph\tremove rejected small contig\n"
        ));
        assert!(apply_text(&text, &mut g, &inputs, &evidence)
            .unwrap_err()
            .to_string()
            .contains("missing node"));
    }
    fn cycle_contract(nodes: &[(&str, &str)], edges: &[(&str, &str)], walk: &[&str])
        -> ResolveContract
    {
        let edge = |a: &str, b: &str| EvidenceLink {
            from: a.to_string(),
            from_orientation: '+',
            to: b.to_string(),
            to_orientation: '+',
            overlap: "0M".to_string(),
            tags: Vec::new(),
        };
        let links = edges.iter().map(|(a, b)| edge(a, b)).collect::<Vec<_>>();
        ResolveContract {
            source_sequences: nodes
                .iter()
                .map(|(n, s)| (n.to_string(), s.to_string()))
                .collect(),
            topology_links: links.clone(),
            resolved_links: links,
            restored_links: Vec::new(),
            lineage: BTreeMap::new(),
            walk: walk
                .iter()
                .map(|n| WalkStep {
                    source_node: n.to_string(),
                    instance_node: n.to_string(),
                    orientation: '+',
                })
                .collect(),
            circular: true,
            circular_evidence: "recorded".to_string(),
            reference_placement_fallback: false,
            alignment_reverse_complemented: false,
            alignment_rotation_step: 0,
            path_issues: Vec::new(),
        }
    }
    #[test]
    fn copy_suffixed_base_node_is_not_treated_as_a_repeat_alias() {
        // Resolve may itself name a merged node "r_copy1"; its duplicate instance is
        // "r_copy1_copy1". Membership-aware resolution must keep "r_copy1" as the base
        // node rather than stripping the suffix to the nonexistent node "r".
        let mut g = cycle_contract(
            &[("a", "AAA"), ("b", "CCC"), ("r_copy1", "GGG")],
            &[("a", "b"), ("b", "r_copy1"), ("r_copy1", "a")],
            &["a", "b", "r_copy1"],
        );
        let (text, inputs, evidence) = record(
            "walk\tpath\ta+,b+,r_copy1+\tr_copy1+,a+,b+\tnative_graph\trotate a circular path that keeps the copy-suffixed base\n",
        );
        let report = apply_text(&text, &mut g, &inputs, &evidence).unwrap();
        assert_eq!(
            g.walk.iter().map(|w| w.source_node.clone()).collect::<Vec<_>>(),
            vec!["r_copy1".to_string(), "a".to_string(), "b".to_string()]
        );
        assert_eq!(report.new_walk, "r_copy1+,a+,b+");
    }
    #[test]
    fn split_node_replaces_linear_parent_with_hash_bound_children_and_reroutes_walk() {
        let mut g = cycle_contract(
            &[("a", "AAAA"), ("p", "ACGTTTGG"), ("b", "GGCC")],
            &[("a", "p"), ("p", "b"), ("b", "a")],
            &["a", "p", "b"],
        );
        let c1 = "ACGT";
        let c2 = "TTGG";
        let fasta = format!(">c1\n{c1}\n>c2\n{c2}\n");
        let path = std::env::temp_dir().join(format!(
            "orgraft_split_test_{}_{}.fasta",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::write(&path, fasta.as_bytes()).unwrap();
        let ev_hash = sha256_hex(fasta.as_bytes());
        let input_hash = |s: &str| sha256_hex(s.as_bytes());
        let mut text = "kind\ttarget\tbefore\tafter\tevidence\tnote\n".to_string();
        for role in ["checked_gfa", "edited_gfa", "resolve_details", "polished_fasta"] {
            text.push_str(&format!(
                "input\t{role}\t{}\t.\t.\tunchanged source\n",
                input_hash(role)
            ));
        }
        text.push_str(&format!(
            "evidence\tcorr\t{ev_hash}\t{}\t.\thash-bound corrected child sequences\n",
            path.display()
        ));
        let parent_before = sha256_hex(b"ACGTTTGG");
        let h1 = sha256_hex(c1.as_bytes());
        let h2 = sha256_hex(c2.as_bytes());
        text.push_str(&format!(
            "split_node\tp>c1,c2\t{parent_before}\t{h1},{h2}\tcorr\tsplit a mis-merged linear parent into two children\n"
        ));
        for (a, b) in [("a", "p"), ("p", "b")] {
            text.push_str(&format!(
                "link\t{a}+>{b}+\t0M\t.\tcorr\tremove the parent's old incident edge\n"
            ));
        }
        for (a, b) in [("a", "c1"), ("c1", "c2"), ("c2", "b")] {
            text.push_str(&format!(
                "link\t{a}+>{b}+\t.\t0M\tcorr\treconnect through the split children\n"
            ));
        }
        text.push_str(
            "walk\tpath\ta+,p+,b+\ta+,c1+,c2+,b+\tcorr\treroute around the split children\n",
        );
        let inputs = ["checked_gfa", "edited_gfa", "resolve_details", "polished_fasta"]
            .into_iter()
            .map(|s| (s.to_string(), input_hash(s)))
            .collect::<BTreeMap<_, _>>();
        let evidence = BTreeMap::from([("corr".to_string(), ev_hash)]);
        let report = apply_text(&text, &mut g, &inputs, &evidence).unwrap();
        assert!(!g.source_sequences.contains_key("p"));
        assert_eq!(g.source_sequences["c1"], c1);
        assert_eq!(g.source_sequences["c2"], c2);
        assert_eq!(report.split_nodes, 1);
        assert_eq!(
            g.walk.iter().map(|w| w.source_node.clone()).collect::<Vec<_>>(),
            vec!["a".to_string(), "c1".to_string(), "c2".to_string(), "b".to_string()]
        );
        for pair in g.walk.windows(2) {
            g.edge_between(&pair[0], &pair[1]).unwrap();
        }
        g.edge_between(g.walk.last().unwrap(), &g.walk[0]).unwrap();
        let _ = fs::remove_file(&path);
    }
}

#[cfg(test)]
mod tests_repeat_complete_graph_span {
    use super::*;

    fn segment(name: &str, sequence: String) -> Segment {
        Segment { name: name.into(), sequence, tags: Vec::new() }
    }
    fn link(a: &str, sa: char, b: &str, sb: char, overlap: usize) -> Link {
        Link { from_name: a.into(), from_orient: sa, to_name: b.into(), to_orient: sb,
            overlap: format!("{overlap}M"), tags: Vec::new() }
    }
    // Dimensions match the affected graph template. The bases are explicit
    // unit fixtures and do not establish an actual sample measurement.
    fn fixture(left_bp: usize, right_bp: usize, left_overlap: usize, right_overlap: usize) -> Gfa {
        let left = format!("{}{}", "A".repeat(left_bp-left_overlap), "G".repeat(left_overlap));
        let right = format!("{}{}", "G".repeat(right_overlap), "T".repeat(right_bp-right_overlap));
        let mut segments = BTreeMap::new();
        for (name, dna) in [("L1", left.clone()), ("L2", left), ("R", "G".repeat(891)),
            ("T1", right.clone()), ("T2", right)] {
            segments.insert(name.into(), segment(name, dna));
        }
        Gfa { headers: Vec::new(), segments,
            order: ["L1","L2","R","T1","T2"].iter().map(|v|v.to_string()).collect(),
            links: vec![link("L1",'+',"R",'+',left_overlap),link("L2",'+',"R",'+',left_overlap),
                link("R",'+',"T1",'+',right_overlap),link("R",'+',"T2",'+',right_overlap)],
            other_lines: Vec::new() }
    }


    fn chain_classes(graph: &Gfa) -> BTreeMap<String,String> {
        graph.segments.keys().map(|n|(n.clone(),if n=="R"{"repeat_node"}else{"single_copy_node"}.into())).collect()
    }
    fn chain_graph() -> Gfa {
        let mut g=fixture(5000,5000,0,0);
        g.segments.get_mut("T1").unwrap().sequence="T".into();
        g.segments.insert("tail".into(),segment("tail","T".repeat(4999)));
        g.order.push("tail".into());g.links.push(link("T1",'+',"tail",'+',0));g
    }
    #[test]
    fn linear_chain_long_arms_and_split_short_arm_preserve_full_template_windows() {
        let original=fixture(5000,5000,0,0);
        let expected=repeat_node_path_templates(&original,"R").unwrap();
        for graph in [&original,&chain_graph()] {
            let actual=repeat_node_linear_chain_templates(graph,&chain_classes(graph),"R",0,0).unwrap();
            assert_eq!(actual.iter().map(|(a,b,d,s,e,_)|(*a,*b,d.clone(),*s,*e)).collect::<Vec<_>>(),expected);
            assert!(actual.iter().all(|(_,_,s,a,b,_)|s.len()==4891&&*a==1000&&*b==3891));
        }
        let graph=chain_graph();
        let old=repeat_node_path_templates(&graph,"R").unwrap();
        assert!(old.iter().any(|(_,_,dna,_,end)|dna.len()==2892&&*end==dna.len()));
        let classes=chain_classes(&graph);
        let a=repeat_path_support_rows_with_linear_chain_flanks(&graph,&classes,&BTreeMap::new(),false,None,Path::new("unused"),1,Path::new("unused"),false).unwrap();
        let b=repeat_path_support_rows(&graph,&classes,None,Path::new("unused"),1,Path::new("unused"),false).unwrap();
        assert_eq!(format!("{a:?}"),format!("{b:?}"));
    }
    #[test]
    fn linear_chain_rejects_branch_cycle_other_repeat_unknown_overlap_and_missing_source() {
        for mode in 0..7 {
            let mut graph=chain_graph();let mut classes=chain_classes(&graph);
            match mode {
                0=>graph.links.push(link("T1",'+',"T2",'+',0)),
                1=>{graph.links.pop();graph.links.push(link("T1",'+',"T1",'+',0));},
                2=>{classes.insert("tail".into(),"repeat_node".into());},
                3=>{classes.remove("tail");},
                4=>graph.links.last_mut().unwrap().overlap="1M".into(),
                5=>{graph.segments.remove("tail");},
                _=>graph.links.push(link("L2",'+',"tail",'+',0)),
            }
            assert!(repeat_node_linear_chain_templates(&graph,&classes,"R",0,0).is_err(),"mode {mode}");
        }
    }
    #[test]
    fn linear_chain_replays_reverse_stored_nodes_and_reciprocal_links() {
        let mut graph=chain_graph();let classes=chain_classes(&graph);
        let expected=repeat_node_linear_chain_templates(&graph,&classes,"R",0,0).unwrap();
        graph.segments.get_mut("tail").unwrap().sequence="A".repeat(4999);
        graph.links.last_mut().unwrap().to_orient='-';
        for edge in &mut graph.links {
            let (a,sa,b,sb)=(edge.from_name.clone(),edge.from_orient,edge.to_name.clone(),edge.to_orient);
            edge.from_name=b;edge.from_orient=flip_orient(sb);edge.to_name=a;edge.to_orient=flip_orient(sa);
        }
        let got=repeat_node_linear_chain_templates(&graph,&classes,"R",0,0).unwrap();
        assert_eq!(got.iter().map(|(a,b,d,s,e,_)|(*a,*b,d,s,e)).collect::<Vec<_>>(),expected.iter().map(|(a,b,d,s,e,_)|(*a,*b,d,s,e)).collect::<Vec<_>>());
    }

    #[test]
    fn complete_span_891_core_and_398_1278_extents_keep_outer_anchor_bounds() {
        let graph=fixture(5000,5000,0,0);
        let templates=repeat_node_path_templates_with_extents(&graph,"R",398,1278).unwrap();
        assert_eq!(templates.len(),4);
        for (_,_,word,start,end) in templates {
            assert_eq!(word,format!("{}{}{}","A".repeat(2398),"G".repeat(891),"T".repeat(3278)));
            assert_eq!(word.len(),6567);
            assert_eq!((start,end),(1000,5567));
            assert_eq!(end-start,891+398+1278+2*REPEAT_PATH_REQUIRED_FLANK_BP);
        }
    }

    #[test]
    fn zero_extents_preserve_default_template_bytes_required_bounds_and_rows() {
        for graph in [fixture(5000,5000,0,0),fixture(300,800,0,0),fixture(5000,5000,17,23)] {
            let old=repeat_node_path_templates(&graph,"R").unwrap();
            assert_eq!(old,repeat_node_path_templates_with_extents(&graph,"R",0,0).unwrap());
            for (a,b,word,start,end) in &old {
                let left=repeat_path_edge(&graph.links[a-1],"R",*a).unwrap();
                let right=repeat_path_edge(&graph.links[b-1],"R",*b).unwrap();
                assert_eq!((word.clone(),*start,*end),spell_repeat_path_template(&graph,"R",&left,&right).unwrap());
            }
            let classes=BTreeMap::from([("R".into(),"repeat_node".into())]);
            let old_rows=repeat_path_support_rows(&graph,&classes,None,Path::new("unused_minimap2"),1,Path::new("unused_debug"),false).unwrap();
            let new_rows=repeat_path_support_rows_with_extents(&graph,&classes,&BTreeMap::new(),None,Path::new("unused_minimap2"),1,Path::new("unused_debug"),false).unwrap();
            assert_eq!(format!("{old_rows:?}"),format!("{new_rows:?}"));
        }
        let short=repeat_node_path_templates(&fixture(300,800,0,0),"R").unwrap();
        assert!(short.iter().all(|(_,_,word,start,end)|word.len()==1991&&*start==0&&*end==1991));
    }

    #[test]
    fn complete_span_deduplicates_nonzero_overlaps_and_both_edge_storage_directions() {
        let graph=fixture(5000,5000,17,23);
        let actual=repeat_node_path_templates_with_extents(&graph,"R",398,1278).unwrap();
        let zero=repeat_node_path_templates_with_extents(&fixture(5000,5000,0,0),"R",398,1278).unwrap();
        assert_eq!(actual,zero);
        let mut reciprocal=graph.clone();
        reciprocal.links=graph.links.iter().map(|edge|Link { from_name:edge.to_name.clone(),
            from_orient:flip_orient(edge.to_orient),to_name:edge.from_name.clone(),
            to_orient:flip_orient(edge.from_orient),overlap:edge.overlap.clone(),tags:edge.tags.clone() }).collect();
        assert_eq!(actual,repeat_node_path_templates_with_extents(&reciprocal,"R",398,1278).unwrap());
        let mut reversed_arms=graph.clone();
        for (name,node) in &mut reversed_arms.segments {
            if name!="R" { node.sequence=reverse_complement(&node.sequence); }
        }
        for edge in &mut reversed_arms.links {
            if edge.from_name!="R" {edge.from_orient=flip_orient(edge.from_orient);}
            if edge.to_name!="R" {edge.to_orient=flip_orient(edge.to_orient);}
        }
        assert_eq!(actual,repeat_node_path_templates_with_extents(&reversed_arms,"R",398,1278).unwrap());
    }

    #[test]
    fn complete_span_requires_every_full_outer_anchor_without_shortening_required_span() {
        let exact=fixture(398+1000,1278+1000,0,0);
        for (_,_,word,start,end) in repeat_node_path_templates_with_extents(&exact,"R",398,1278).unwrap() {
            assert_eq!((word.len(),start,end),(4567,0,4567));
        }
        assert!(repeat_node_path_templates_with_extents(&fixture(1397,2278,0,0),"R",398,1278).is_err());
        assert!(repeat_node_path_templates_with_extents(&fixture(1398,2277,0,0),"R",398,1278).is_err());
        let classes=BTreeMap::from([("R".into(),"repeat_node".into())]);
        let extents=BTreeMap::from([("R".into(),(398,1278))]);
        assert!(repeat_path_support_rows_with_extents(&fixture(1397,2278,0,0),&classes,&extents,None,
            Path::new("unused_minimap2"),1,Path::new("unused_debug"),false).is_err());
    }

    #[test]
    fn complete_span_rejects_extent_and_overlap_overflow_before_template_slicing() {
        let graph=fixture(5000,5000,0,0);
        assert!(repeat_node_path_templates_with_extents(&graph,"R",usize::MAX,1278).is_err());
        assert!(repeat_node_path_templates_with_extents(&graph,"R",398,usize::MAX).is_err());
        let mut invalid=graph.clone();invalid.links[0].overlap=format!("{}M1M",usize::MAX);
        assert!(repeat_node_path_templates_with_extents(&invalid,"R",398,1278).is_err());
        invalid.links[0].overlap="17X".into();
        assert!(repeat_node_path_templates_with_extents(&invalid,"R",398,1278).is_err());
    }

    #[test]
    fn complete_span_rejects_wrong_overlap_dna_missing_nodes_and_unbound_extent_keys() {
        let mut graph=fixture(5000,5000,17,23);graph.segments.get_mut("L1").unwrap().sequence.replace_range(4999..5000,"A");
        assert!(repeat_node_path_templates_with_extents(&graph,"R",398,1278).is_err());
        graph.links[0].overlap="892M".into();
        assert!(repeat_node_path_templates_with_extents(&graph,"R",398,1278).is_err());
        let graph=fixture(5000,5000,0,0);let classes=BTreeMap::from([("R".into(),"repeat_node".into())]);
        let mut no_left=graph.clone();no_left.links.retain(|edge|edge.from_name=="R");
        let extents=BTreeMap::from([("R".into(),(398,1278))]);
        assert!(repeat_path_support_rows_with_extents(&no_left,&classes,&extents,None,
            Path::new("unused_minimap2"),1,Path::new("unused_debug"),false).is_err());
        for name in ["absent","L1"] {
            let extents=BTreeMap::from([(name.into(),(398,1278))]);
            assert!(repeat_path_support_rows_with_extents(&graph,&classes,&extents,None,
                Path::new("unused_minimap2"),1,Path::new("unused_debug"),false).is_err());
        }
    }
}

#[cfg(test)]
mod tests_restored_projection {
    use super::*;

    fn unique_dna(length: usize) -> String {
        let mut state = 37u64;
        (0..length).map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state >> 62) as usize] as char
        }).collect()
    }

    #[test]
    fn gap_rich_single_copy_keeps_original_substitution_and_mapping_gates() {
        // The saved projection regression's complete PAF spans and EQX counts.
        let line = "arm\t86458\t0\t86458\t+\tcircle\t400000\t82977\t166729\t83747\t86460\t60\ttp:A:P\tcg:Z:40000=2708I2D3X43747=";
        let hit = parse_paf(line).unwrap();
        assert!(hit.identity() < 0.97);
        assert!(usable_node_projection_hit(&hit, false));
        assert!(!usable_node_projection_hit(&hit, true));
        for mapq in [0, 19] {
            let mut low_mapq = hit.clone();
            low_mapq.mapq = mapq;
            assert!(!usable_node_projection_hit(&low_mapq, false));
        }
        let mut short = hit.clone();
        short.query_end = 86458 * 3 / 4 - 1;
        assert!(!usable_node_projection_hit(&short, false));
        let mut substitutions = hit.clone();
        substitutions.substitution_mismatches = 1000;
        assert!(!usable_node_projection_hit(&substitutions, false));
    }

    #[test]
    fn doubled_circle_requires_exact_unique_full_query_pair() {
        let line = |start, strand, primary, cigar: &str| format!(
            "arm\t250\t0\t250\t{strand}\tcircle\t2000\t{start}\t{}\t250\t270\t0\ttp:A:{primary}\tcg:Z:{cigar}", start + 270);
        let first = parse_paf(&line(100, '-', 'P', "120=20D130=")).unwrap();
        let second = parse_paf(&line(1100, '-', 'S', "120=20D130=")).unwrap();
        assert!(!usable_node_projection_hit(&first, false));
        assert_eq!(unique_circular_gap_projection_hit(&[first.clone(), second.clone()], 1000).unwrap().target_start, 100);
        let mut variants = Vec::new();
        let mut x = second.clone(); x.strand = '+'; variants.push(x);
        let mut x = second.clone(); x.primary = true; variants.push(x);
        let mut x = second.clone(); x.target_start += 1; variants.push(x);
        let mut x = second.clone(); x.target_end -= 1; variants.push(x);
        let mut x = second.clone(); x.query_end -= 1; variants.push(x);
        let mut x = second.clone(); x.mapq = 1; variants.push(x);
        let mut x = second.clone(); x.substitution_mismatches = 10; variants.push(x);
        variants.push(parse_paf(&line(1100, '-', 'S', "120M20D130M")).unwrap());
        variants.push(parse_paf(&line(1100, '-', 'S', "99=20D151=")).unwrap());
        for invalid in variants {
            assert!(unique_circular_gap_projection_hit(&[first.clone(), invalid], 1000).is_none());
        }
        assert!(unique_circular_gap_projection_hit(&[first.clone()], 1000).is_none());
        assert!(unique_circular_gap_projection_hit(&[first.clone(), second.clone(), second], 1000).is_none());
    }

    #[test]
    fn unsupported_single_copy_does_not_claim_another_nodes_gap() {
        let graph = Gfa { headers: vec![], segments: BTreeMap::new(),
            order: vec!["unsupported".into(), "supported".into()],
            links: vec![], other_lines: vec![] };
        let classes = graph.order.iter().map(|x| (x.clone(), "single_copy_node".into())).collect();
        let hit = make_synthetic_hit("supported", 100, 1, 100, 100, 199, '+');
        let mapping = Mapping {
            by_node: BTreeMap::from([
                ("unsupported".into(), NodeMapping { selected: None, accepted_hits: vec![], all_hits: vec![] }),
                ("supported".into(), NodeMapping { selected: Some(hit.clone()), accepted_hits: vec![hit.clone()], all_hits: vec![hit] }),
            ]), reference_len: 1000, command: String::new(), stderr: String::new(),
        };
        let gaps = vec![RepeatGap { left_repeat_node: "r".into(), left_copy_index: 1, right_repeat_node: "r".into(), right_copy_index: 2, parts: vec![(100, 199)] }];
        let assigned = assign_single_copy_nodes_to_gaps(&graph, &classes, &mapping, &gaps);
        assert!(!assigned.contains_key("unsupported"));
        assert_eq!(assigned["supported"].parts, vec![(100, 199)]);
    }

    fn seam_fixture() -> (String, String, String) {
        // Real full-query PAF geometry: 287102 query bases, 2232 unaligned
        // circular target bases, and one query-only base in the aligned core.
        let target = unique_dna(289333);
        let frame = rotate(&target, 273503);
        let source = format!("{}A{}", &frame[..282021], &frame[282021..287101]);
        let paf = "path_source\t287102\t0\t287102\t+\tpath_target\t578666\t273503\t560604\t287101\t287102\t60\ttp:A:P\tcg:Z:282021=1I5080=".to_string();
        (source, target, paf)
    }

    #[test]
    fn saved_complete_seam_geometry_replays_all_target_bases_and_graph() {
        let (source, target, paf) = seam_fixture();
        let candidates = candidates_from_paf(&source, &target, true, &paf).unwrap();
        assert_eq!(candidates.len(), 1);
        let candidate = &candidates[0];
        assert_eq!(candidate.method, "minimap2_primary_complete_path_circular_seam_gap");
        assert_eq!(candidate.cigar, "282021=1I5080=2232D");
        assert_eq!(candidate.target_in_path_frame.len(), target.len());
        let nodes = BTreeMap::from([("node".to_string(), source)]);
        let edge = PathEdge { from: "node".into(), from_reverse: false, to: "node".into(), to_reverse: false, overlap: 0 };
        let walk = vec![PathStep { node: "node".into(), reverse: false }];
        let prepared = prepare_path(&nodes, &[edge.clone()], &walk, true).unwrap();
        let projected = project_path_sequence(&prepared, &candidate.target_in_path_frame, &candidate.cigar, 100).unwrap();
        let replay = prepare_path(&projected.sequences, &[edge], &walk, true).unwrap();
        assert_eq!(replay.source_spell, candidate.target_in_path_frame);
    }

    #[test]
    fn seam_candidate_refuses_wrong_or_nonunique_anchor_evidence() {
        let (source, target, paf) = seam_fixture();
        let anchor = parse_anchor(&paf, source.len(), target.len() * 2).unwrap();
        assert!(circular_seam_gap_candidate(&source, &target, &anchor, 2232).is_some());
        for changed in [paf.replace("\t+\tpath_target", "\t-\tpath_target"), paf.replace("cg:Z:282021=1I5080=", "cg:Z:282020=1I5080=")] {
            assert!(parse_anchor(&changed, source.len(), target.len() * 2).is_err());
        }
        for changed in [paf.replace("\t60\t", "\t19\t"), paf.replace("tp:A:P", "tp:A:S"), paf.replace("282021=1I5080=", "282021M1I5080M")] {
            let invalid = parse_anchor(&changed, source.len(), target.len() * 2).unwrap();
            assert!(circular_seam_gap_candidate(&source, &target, &invalid, 2232).is_none());
        }
        assert!(circular_seam_gap_candidate(&source, &target, &anchor, 0).is_none());
        assert!(circular_seam_gap_candidate(&source, &target, &anchor, MAX_TERMINAL_TARGET_BP + 1).is_none());
        let repetitive_target = "A".repeat(target.len());
        let repetitive_source = "A".repeat(source.len());
        assert!(circular_seam_gap_candidate(&repetitive_source, &repetitive_target, &anchor, 2232).is_none());
        let mut wrong = source.clone().into_bytes();
        wrong[120] = if wrong[120] == b'A' { b'C' } else { b'A' };
        assert!(circular_seam_gap_candidate(std::str::from_utf8(&wrong).unwrap(), &target, &anchor, 2232).is_none());
        let mut short_query = anchor;
        short_query.query_start = 1;
        assert!(circular_seam_gap_candidate(&source, &target, &short_query, 2232).is_none());
    }
}

#[cfg(test)]
mod tests_circular_seam_overlap {
    use super::*;

    fn unique_dna(n: usize, seed: u64) -> String {
        let mut state = seed;
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize] as char
            })
            .collect()
    }

    fn fixture() -> (PreparedPath, String, String, BTreeMap<String, String>) {
        let seam = unique_dna(251, 0xabcdef);
        let a = format!("{seam}{}{seam}", unique_dna(809, 0x987654));
        let r = format!("{seam}{}{seam}", unique_dna(617, 0x123456));
        let b = format!("{seam}{}{seam}", unique_dna(733, 0xfedcba));
        let nodes = BTreeMap::from([
            ("a".into(), a.clone()),
            ("r".into(), r.clone()),
            ("b".into(), reverse_complement_iupac(&b)),
        ]);
        let walk = vec![
            PathStep {
                node: "a".into(),
                reverse: false,
            },
            PathStep {
                node: "r".into(),
                reverse: false,
            },
            PathStep {
                node: "b".into(),
                reverse: true,
            },
            PathStep {
                node: "r".into(),
                reverse: false,
            },
        ];
        let edges = (0..walk.len())
            .map(|i| {
                let left = &walk[i];
                let right = &walk[(i + 1) % walk.len()];
                PathEdge {
                    from: left.node.clone(),
                    from_reverse: left.reverse,
                    to: right.node.clone(),
                    to_reverse: right.reverse,
                    overlap: 0,
                }
            })
            .collect::<Vec<_>>();
        let prepared = prepare_path(&nodes, &edges, &walk, true).unwrap();
        let target = format!("{a}{}{}{}", &r[251..], &b[251..], &r[251..r.len() - 251]);
        let cigar = format!(
            "{}=251I{}=251I{}=251I{}=",
            a.len(),
            r.len() - 251,
            b.len() - 251,
            r.len() - 251
        );
        let paf = format!(
            "path_source\t{}\t0\t{}\t+\tpath_target\t{}\t0\t{}\t{}\t{}\t60\ttp:A:P\tcg:Z:{cigar}",
            prepared.source_spell.len(),
            prepared.source_spell.len(),
            target.len() * 2,
            target.len() + 251,
            target.len() + 251,
            prepared.source_spell.len()
        );
        (prepared, target, paf, nodes)
    }

    #[test]
    fn complete_unfolded_anchor_preserves_all_copies_edges_and_query() {
        let (prepared, target, paf, nodes) = fixture();
        let source_before = prepared.source_spell.clone();
        let candidates = candidates_from_paf(&prepared.source_spell, &target, true, &paf).unwrap();
        assert_eq!(candidates.len(), 1);
        let candidate = &candidates[0];
        assert_eq!(candidate.terminal_target_overlap, 251);
        let ops = parse_cigar(&candidate.cigar).unwrap();
        assert_eq!(
            consumed(&ops).unwrap(),
            (source_before.len(), target.len() + 251)
        );
        let plan = project_circular_seam_overlap(&prepared, candidate).unwrap();
        assert_eq!(
            plan.sequences.keys().collect::<Vec<_>>(),
            nodes.keys().collect::<Vec<_>>()
        );
        assert_eq!(plan.copies.len(), 4);
        assert_eq!(
            plan.copies.iter().filter(|copy| copy.node == "r").count(),
            2
        );
        assert_eq!(plan.edge_overlaps.len(), prepared.edges.len());
        assert_eq!(
            plan.edge_overlaps
                .iter()
                .filter(|edge| edge.new_overlap == 251)
                .count(),
            1
        );
        let edges = prepared
            .edges
            .iter()
            .cloned()
            .zip(&plan.edge_overlaps)
            .map(|(mut edge, overlap)| {
                edge.overlap = overlap.new_overlap;
                edge
            })
            .collect::<Vec<_>>();
        let walk = prepared
            .occurrences
            .iter()
            .map(|copy| copy.step.clone())
            .collect::<Vec<_>>();
        let replay = prepare_path(&plan.sequences, &edges, &walk, true).unwrap();
        assert_eq!(replay.source_spell, rotate(&target, 251));
        assert_eq!(replay.source_copies(), plan.copies);
        assert_eq!(prepared.source_spell, source_before);
        assert_eq!(
            prepared
                .edges
                .iter()
                .map(|edge| edge.overlap)
                .collect::<Vec<_>>(),
            vec![0; 4]
        );
    }

    #[test]
    fn incomplete_wrong_strand_nonprimary_and_low_quality_witnesses_are_rejected() {
        let (prepared, target, paf, _) = fixture();
        for changed in [
            paf.replace("\t+\tpath_target", "\t-\tpath_target"),
            paf.replace("\t60\t", "\t19\t"),
            paf.replace("tp:A:P", "tp:A:S"),
            paf.replace("1311=", "1311M"),
        ] {
            assert!(candidates_from_paf(&prepared.source_spell, &target, true, &changed).is_err());
        }
        assert!(candidates_from_paf(&prepared.source_spell, &target, false, &paf).is_err());
        let mut anchor = parse_anchor(&paf, prepared.source_spell.len(), target.len() * 2).unwrap();
        anchor.query_start = 1;
        assert!(
            circular_seam_overlap_candidate(&prepared.source_spell, &target, &anchor).is_none()
        );
        anchor.query_start = 0;
        let mut wrong = prepared.source_spell.clone().into_bytes();
        wrong[700] = if wrong[700] == b'A' { b'C' } else { b'A' };
        assert!(
            circular_seam_overlap_candidate(std::str::from_utf8(&wrong).unwrap(), &target, &anchor)
                .is_none()
        );
        let mut wrong_end = prepared.source_spell.clone().into_bytes();
        let last = *wrong_end.last().unwrap();
        *wrong_end.last_mut().unwrap() = if last == b'A' { b'C' } else { b'A' };
        assert!(
            circular_seam_overlap_candidate(
                std::str::from_utf8(&wrong_end).unwrap(),
                &target,
                &anchor
            )
            .is_none()
        );
    }

    #[test]
    fn shared_core_projection_uses_the_complete_unfolded_last_copy() {
        let (prepared, target, paf, _) = fixture();
        let alternate = paf.replace("1311=251I868=", "1060=251I1119=");
        let candidate = candidates_from_paf(&prepared.source_spell, &target, true, &alternate)
            .unwrap()
            .remove(0);
        let plan = project_circular_seam_overlap(&prepared, &candidate).unwrap();
        assert!(!plan.shared_core_trims.is_empty());
        assert_eq!(
            plan.copies.iter().filter(|copy| copy.node == "r").count(),
            2
        );
        let edges = prepared
            .edges
            .iter()
            .cloned()
            .zip(&plan.edge_overlaps)
            .map(|(mut edge, overlap)| {
                edge.overlap = overlap.new_overlap;
                edge
            })
            .collect::<Vec<_>>();
        let walk = prepared
            .occurrences
            .iter()
            .map(|copy| copy.step.clone())
            .collect::<Vec<_>>();
        let replay = prepare_path(&plan.sequences, &edges, &walk, true).unwrap();
        assert_eq!(
            replay.source_spell,
            rotate(
                &target,
                plan.target_frame_shift.rem_euclid(target.len() as isize) as usize
            )
        );
        assert_eq!(replay.source_copies(), plan.copies);
        assert_eq!(
            consumed(&parse_cigar(&plan.alignment_cigar).unwrap()).unwrap(),
            (prepared.source_spell.len(), target.len() + 251)
        );
    }

    #[test]
    fn equally_supported_shared_cores_remain_unresolved() {
        let nodes = BTreeMap::from([("r".into(), "AAAA".into())]);
        let edge = PathEdge {
            from: "r".into(),
            from_reverse: false,
            to: "r".into(),
            to_reverse: false,
            overlap: 0,
        };
        let walk = vec![
            PathStep {
                node: "r".into(),
                reverse: false
            };
            2
        ];
        let prepared = prepare_path(&nodes, &[edge], &walk, true).unwrap();
        let target = b"ACGCATG";
        let domains = cigar_domains(prepared.source_spell.as_bytes(), b"ACGCATGA", "8M").unwrap();
        assert!(
            project_shared_cores_with_terminal_overlap(&prepared, target, &domains, 1)
                .unwrap_err()
                .contains("shared_core_is_not_unique")
        );
    }

    #[test]
    fn repeated_complete_circle_is_not_a_unique_frame() {
        let target = "ACGT".repeat(1000);
        let source = format!("{target}{}", &target[..251]);
        let anchor = Anchor {
            query_start: 0,
            query_end: source.len(),
            target_start: 0,
            target_end: source.len(),
            cigar: vec![Op {
                count: source.len(),
                code: b'=',
            }],
            mapq: 60,
            primary: true,
        };
        assert!(circular_seam_overlap_candidate(&source, &target, &anchor).is_none());
    }

    #[test]
    fn seam_projection_keeps_all_edge_and_occurrence_constraints() {
        let (prepared, target, paf, _) = fixture();
        let candidate = candidates_from_paf(&prepared.source_spell, &target, true, &paf)
            .unwrap()
            .remove(0);
        let mut linear = prepared.clone();
        linear.circular = false;
        assert!(project_circular_seam_overlap(&linear, &candidate).is_err());
        let mut conflicting = prepared.clone();
        conflicting.edges.push(PathEdge {
            from: "a".into(),
            from_reverse: false,
            to: "b".into(),
            to_reverse: true,
            overlap: 252,
        });
        assert!(
            project_circular_seam_overlap(&conflicting, &candidate)
                .unwrap_err()
                .contains("edge_overlap_incompatible")
        );
        let mut changed = candidate.clone();
        changed.cigar.push_str("1=");
        assert!(project_circular_seam_overlap(&prepared, &changed).is_err());
        let mut missing = prepared.clone();
        missing.occurrences.pop();
        assert!(project_circular_seam_overlap(&missing, &candidate).is_err());
    }

    #[test]
    fn two_distinct_complete_witnesses_do_not_choose_a_convenient_boundary() {
        let (prepared, target, paf, _) = fixture();
        let alternative = paf.replace("1311=251I868=", "1060=251I1119=");
        assert_ne!(alternative, paf);
        let anchor =
            parse_anchor(&alternative, prepared.source_spell.len(), target.len() * 2).unwrap();
        assert!(
            circular_seam_overlap_candidate(&prepared.source_spell, &target, &anchor).is_some()
        );
        let error = candidates_from_paf(
            &prepared.source_spell,
            &target,
            true,
            &format!("{paf}\n{alternative}"),
        )
        .unwrap_err();
        assert!(error.contains("multiple_complete_alignment_witnesses"));
    }

    #[test]
    fn a_target_insertion_at_an_original_node_cut_is_not_chosen_arbitrarily() {
        let (prepared, target, paf, _) = fixture();
        let mut candidate = candidates_from_paf(&prepared.source_spell, &target, true, &paf)
            .unwrap()
            .remove(0);
        let cut = 1311;
        candidate.target_in_path_frame.insert(cut, 'A');
        candidate.cigar = candidate.cigar.replacen("1311=", "1311=1D", 1);
        assert!(
            project_circular_seam_overlap(&prepared, &candidate)
                .unwrap_err()
                .contains("ambiguous_node_boundary")
        );
    }
}


// A validated ordinary checkpoint may have changed the selected repeat pairing.
// Reproduce that native transition from its original inputs; never use the mutable
// historical rebuild output or its old VS coordinates as the final authority.
mod ordinary_source {
    use super::*;
    type Values = BTreeMap<String, String>;
    fn err(message: impl Into<String>) -> OrgraftError {
        OrgraftError::InvalidArgument(format!("ordinary source transition: {}", message.into()))
    }
    fn require(ok: bool, message: &str) -> Result<(), OrgraftError> {
        if ok {
            Ok(())
        } else {
            Err(err(message))
        }
    }
    fn fields(path: &Path, width: usize) -> Result<Vec<Vec<String>>, OrgraftError> {
        let text = fs::read_to_string(path)?;
        let mut lines = text.lines();
        require(lines.next().is_some(), "empty evidence table")?;
        lines
            .filter(|l| !l.is_empty())
            .map(|line| {
                let f = line.split('\t').map(str::to_string).collect::<Vec<_>>();
                require(f.len() == width, "malformed evidence row")?;
                Ok(f)
            })
            .collect()
    }
    fn metrics(path: &Path) -> Result<Values, OrgraftError> {
        let mut result = Values::new();
        for f in fields(path, 2)? {
            require(
                result.insert(f[0].clone(), f[1].clone()).is_none(),
                "duplicate evidence key",
            )?;
        }
        Ok(result)
    }
    fn report(path: &Path) -> Result<Values, OrgraftError> {
        let mut result = Values::new();
        for f in fields(path, 4)? {
            require(
                result
                    .insert(format!("{}/{}/{}", f[0], f[1], f[2]), f[3].clone())
                    .is_none(),
                "duplicate native report key",
            )?;
        }
        Ok(result)
    }
    fn get<'a>(values: &'a Values, key: &str) -> Result<&'a str, OrgraftError> {
        values
            .get(key)
            .map(String::as_str)
            .ok_or_else(|| err(format!("missing {key}")))
    }
    fn number(values: &Values, key: &str) -> Result<usize, OrgraftError> {
        get(values, key)?
            .parse()
            .map_err(|_| err(format!("invalid {key}")))
    }
    fn same_path(a: &Path, b: &Path) -> Result<(), OrgraftError> {
        require(
            fs::canonicalize(a)? == fs::canonicalize(b)?,
            "evidence path differs from its declared native input",
        )
    }
    fn hash(path: &Path) -> Result<String, OrgraftError> {
        crate::commands::refine::sha256_reader(&mut File::open(path)?)
    }
    fn same_file(a: &Path, b: &Path) -> Result<(), OrgraftError> {
        if fs::canonicalize(a)? == fs::canonicalize(b)? {
            return Ok(());
        }
        require(hash(a)? == hash(b)?, "native input content differs")
    }
    fn one(path: &Path) -> Result<FastaRecord, OrgraftError> {
        let mut records = read_fasta(path)?;
        require(
            records.len() == 1 && !records[0].sequence.is_empty(),
            "expected one nonempty native FASTA record",
        )?;
        Ok(records.remove(0))
    }
    fn equal_fasta(a: &FastaRecord, b: &FastaRecord) -> Result<(), OrgraftError> {
        require(
            a.id == b.id && a.sequence.eq_ignore_ascii_case(&b.sequence),
            "complete FASTA identity/DNA does not replay",
        )
    }
    fn bind(files: &mut BTreeMap<PathBuf, String>, path: &Path) -> Result<(), OrgraftError> {
        let p = fs::canonicalize(path)?;
        let h = hash(&p)?;
        if let Some(old) = files.insert(p, h.clone()) {
            require(old == h, "source changed during replay")?;
        }
        Ok(())
    }
    fn stable(files: &BTreeMap<PathBuf, String>) -> Result<(), OrgraftError> {
        for (p, h) in files {
            require(hash(p)? == *h, "bound source changed")?;
        }
        Ok(())
    }
    pub(super) struct Context {
        round: PathBuf,
        checkpoint: PathBuf,
        source: PathBuf,
        first_report: PathBuf,
        repeat: String,
        source_round: usize,
        final_round: usize,
        files: BTreeMap<PathBuf, String>,
    }
    pub(crate) fn checkpoint(checked: &Path) -> Result<PathBuf, OrgraftError> {
        let checked = fs::canonicalize(checked)?;
        require(
            checked.file_name().and_then(|n| n.to_str()) == Some("checked_draft.gfa")
                && checked
                    .parent()
                    .and_then(Path::file_name)
                    .and_then(|n| n.to_str())
                    == Some("checkpoint_1"),
            "ordinary transition requires the native checked checkpoint locator",
        )?;
        Ok(checked
            .parent()
            .unwrap()
            .parent()
            .ok_or_else(|| err("missing checkpoint owner"))?
            .join("checkpoint_2"))
    }
    fn context(
        round: &Path,
        checked: &Path,
        component: &str,
        organelle: &str,
        reads: &Path,
        details: &Path,
    ) -> Result<Context, OrgraftError> {
        let round = fs::canonicalize(round)?;
        let final_round = round
            .file_name()
            .and_then(|n| n.to_str())
            .and_then(|s| s.strip_prefix("round_"))
            .and_then(|s| s.parse::<usize>().ok())
            .filter(|n| (2..=10).contains(n))
            .ok_or_else(|| err("invalid final ordinary round"))?;
        let owner = round
            .parent()
            .ok_or_else(|| err("missing ordinary owner"))?;
        require(
            owner.file_name().and_then(|n| n.to_str()) == Some(component)
                && owner
                    .parent()
                    .and_then(Path::file_name)
                    .and_then(|n| n.to_str())
                    == Some(organelle),
            "ordinary source identity differs",
        )?;
        let checkpoint = checkpoint(checked)?;
        require(
            checkpoint
                .parent()
                .and_then(Path::file_name)
                .and_then(|n| n.to_str())
                == Some(component)
                && checkpoint
                    .parent()
                    .and_then(Path::parent)
                    .and_then(Path::file_name)
                    .and_then(|n| n.to_str())
                    == Some(organelle),
            "checkpoint identity differs",
        )?;
        let marker_path = checkpoint.join("repeat_pairing.status.tsv");
        let marker = metrics(&marker_path)?;
        let source_round = number(&marker, "source_round")?;
        // The initial implementation handles one native transition followed by ordinary edits.
        require(
            source_round == 1
                && number(&marker, "validation_round")? == final_round
                && get(&marker, "state")? == "validated",
            "requires one validated first-round repeat transition and its final ordinary round",
        )?;
        let repeat = get(&marker, "repeat_node")?.to_string();
        let source_dir = owner.join(format!("round_{source_round}"));
        let first_report = source_dir.join("logs/report.tsv");
        let first = report(&first_report)?;
        let source = source_dir.join("02.polish/polished_aln.fasta");
        same_path(
            Path::new(get(&first, "output/validation_fasta/path")?),
            &source,
        )?;
        let initial =
            metrics(&checkpoint.join(format!("round_{source_round}/checkpoint_2.status.tsv")))?;
        require(
            get(&initial, "status")? == "rebuild_ready"
                && get(&initial, "sv_correction_mode")? == "repeat_pairing_rebuild"
                && get(&initial, "repeat_node")? == repeat
                && number(&initial, "repeat_pairing_source_round")? == source_round
                && number(&initial, "repeat_pairing_validation_round")? == source_round + 1,
            "original checkpoint did not request this native repeat transition",
        )?;
        same_path(
            Path::new(get(&initial, "repeat_pairing_marker")?),
            &marker_path,
        )?;
        let mismatch = checkpoint.join(format!(
            "round_{source_round}/repeat_pairing_repair/repeat_pairing_mismatch.tsv"
        ));
        same_path(Path::new(get(&marker, "mismatch_report")?), &mismatch)?;
        same_path(
            Path::new(get(&initial, "repeat_pairing_mismatch_report")?),
            &mismatch,
        )?;
        // The mutable old rebuild is a historical locator only, never read as sequence.
        require(
            get(&initial, "corrected_fasta")? == get(&marker, "rebuild_fasta")?,
            "historical rebuild locator changed",
        )?;
        let next = report(&owner.join(format!("round_{}/logs/report.tsv", source_round + 1)))?;
        require(
            get(&next, "input/draft/path")? == get(&marker, "rebuild_fasta")?,
            "validation did not consume the historical rebuild transition",
        )?;
        let mismatch_rows = metrics(&mismatch)?;
        require(
            get(&mismatch_rows, "status")? == "repeat_pairing_mismatch"
                && get(&mismatch_rows, "repeat_node")? == repeat,
            "repeat trigger is absent or ambiguous",
        )?;
        let mut files = BTreeMap::new();
        for p in [
            checked,
            details,
            reads,
            &marker_path,
            &first_report,
            &source,
            &mismatch,
            &checkpoint.join(format!("round_{source_round}/checkpoint_2.status.tsv")),
        ] {
            bind(&mut files, p)?;
        }
        for n in source_round..=final_round {
            let dir = owner.join(format!("round_{n}"));
            let rp = dir.join("logs/report.tsv");
            let r = report(&rp)?;
            for (key, value) in [
                ("global/run/organelle", organelle.to_string()),
                ("global/run/subgraph", component.to_string()),
                ("global/run/validate_round", n.to_string()),
                ("global/run/per_read_variant_calls", "true".into()),
            ] {
                require(
                    get(&r, key)? == value,
                    "native report identity or read-level validation differs",
                )?;
            }
            same_path(Path::new(get(&r, "output/workflow_round/dir")?), &dir)?;
            for stage in ["sv-eval", "snv-indel-eval"] {
                require(
                    get(&r, &format!("stage/{stage}:round_{n}/status"))? == "ok",
                    "native validation did not complete",
                )?;
            }
            let fa = if n == source_round {
                source.clone()
            } else {
                dir.join(format!("01.inputs/linear_subgraph.round_{n}.fasta"))
            };
            same_path(Path::new(get(&r, "output/validation_fasta/path")?), &fa)?;
            require(
                one(&fa)?.id == component,
                "ordinary FASTA record ID differs",
            )?;
            same_file(Path::new(get(&r, "input/reads/path")?), reads)?;
            same_file(&dir.join("01.inputs/subgraph_reads.fastq.gz"), reads)?;
            same_file(
                Path::new(get(&r, "input/reference/path")?),
                Path::new(get(&first, "input/reference/path")?),
            )?;
            same_file(
                &dir.join("01.inputs/rotated_reference.fasta"),
                Path::new(get(&first, "input/reference/path")?),
            )?;
            for p in [
                &rp,
                &fa,
                Path::new(get(&r, "input/reads/path")?),
                Path::new(get(&r, "input/reference/path")?),
                &dir.join("01.inputs/subgraph_reads.fastq.gz"),
                &dir.join("01.inputs/rotated_reference.fasta"),
            ] {
                bind(&mut files, p)?;
            }
            if n > source_round {
                require(
                    get(&r, "sv_eval/round_1/sv_support_status")? == "pass",
                    "ordinary validation SV gate failed",
                )?;
                let summary = dir.join("03.validate/03.reports/sv_snv_indel_summary.tsv");
                require(
                    get(&metrics(&summary)?, "status")? == "pass",
                    "ordinary SV summary failed",
                )?;
                bind(&mut files, &summary)?;
            }
        }
        // The original checked graph and ordinary input are independently tied by resolve provenance.
        let ds = fs::read_to_string(details)?;
        let mut resolved = None;
        let mut checked_record = None;
        for l in ds.lines().skip(1) {
            let f = l.split('\t').collect::<Vec<_>>();
            if f.len() < 4 {
                continue;
            }
            if f[0] == "run" && f[2] == "checked_draft_gfa" {
                checked_record = Some(f[3].to_string());
            }
            if f[0] == "run" && f[2] == "resolved_fasta" {
                resolved = Some(f[3].to_string());
            }
        }
        same_file(
            Path::new(
                checked_record
                    .as_deref()
                    .ok_or_else(|| err("resolve lacks original checked input"))?,
            ),
            checked,
        )?;
        let draft = PathBuf::from(get(&first, "input/draft/path")?);
        same_path(
            &draft,
            Path::new(
                resolved
                    .as_deref()
                    .ok_or_else(|| err("resolve lacks original draft"))?,
            ),
        )?;
        bind(&mut files, &draft)?;
        bind(&mut files, Path::new(checked_record.as_ref().unwrap()))?;
        Ok(Context {
            round,
            checkpoint,
            source,
            first_report,
            repeat,
            source_round,
            final_round,
            files,
        })
    }
    fn replay_edits(c: &mut Context, initial: &FastaRecord) -> Result<(), OrgraftError> {
        let owner = c.round.parent().unwrap();
        let mut previous = one(&owner.join(format!(
            "round_{}/01.inputs/linear_subgraph.round_{}.fasta",
            c.source_round + 1,
            c.source_round + 1
        )))?;
        equal_fasta(initial, &previous)?;
        for n in c.source_round + 1..=c.final_round {
            let dir = owner.join(format!("round_{n}"));
            let status_path = c
                .checkpoint
                .join(format!("round_{n}/checkpoint_2.status.tsv"));
            let status = metrics(&status_path)?;
            let high = dir.join("03.validate/03.reports/snv_indel_high.tsv");
            same_path(Path::new(get(&status, "snv_indel_high")?), &high)?;
            same_path(
                Path::new(get(&status, "sv_summary")?),
                &dir.join("03.validate/03.reports/sv_snv_indel_summary.tsv"),
            )?;
            bind(&mut c.files, &status_path)?;
            bind(&mut c.files, &high)?;
            if n == c.final_round {
                require(
                    get(&status, "status")? == "complete"
                        && crate::workflow::ordinary_validation_edit_count(&high)? == 0,
                    "final ordinary checkpoint is not complete",
                )?;
                break;
            }
            require(
                get(&status, "status")? == "next_round_ready"
                    && get(&status, "correction_kind")? == "snv_indel",
                "intermediate checkpoint is not an ordinary edit",
            )?;
            let edits = c.checkpoint.join(format!("round_{n}/pos_ref_alt.txt"));
            same_path(Path::new(get(&status, "pos_ref_alt")?), &edits)?;
            let corrected =
                crate::workflow::replay_ordinary_source_edits(&high, &edits, &previous.sequence)?;
            let next_dir = owner.join(format!("round_{}", n + 1));
            let next =
                one(&next_dir.join(format!("01.inputs/linear_subgraph.round_{}.fasta", n + 1)))?;
            let r = report(&next_dir.join("logs/report.tsv"))?;
            require(
                get(&r, "input/draft/path")? == get(&status, "corrected_fasta")?,
                "ordinary validation input does not follow its checkpoint",
            )?;
            equal_fasta(
                &FastaRecord {
                    header: previous.header.clone(),
                    id: previous.id.clone(),
                    sequence: corrected,
                },
                &next,
            )?;
            bind(&mut c.files, &edits)?;
            previous = next;
        }
        Ok(())
    }
    fn support(
        path: &Path,
        original: &Gfa,
        classes: &BTreeMap<String, String>,
        root: &Path,
    ) -> Result<Vec<RepeatPathSupportRow>, OrgraftError> {
        let mut expected =
            repeat_path_support_rows(original, classes, None, Path::new("unused"), 1, root, false)?;
        let rows = fields(path, 8)?;
        require(
            rows.len() == expected.len(),
            "native support table omitted or added a path",
        )?;
        let mut per_repeat = BTreeMap::<String, BTreeSet<String>>::new();
        let mut counters = BTreeMap::<String, usize>::new();
        for (f, row) in rows.iter().zip(&mut expected) {
            let i = counters.entry(row.repeat_node.clone()).or_default();
            *i += 1;
            require(
                f[0] == row.repeat_node
                    && f[1] == format!("p{i}")
                    && f[2] == row.left_endpoint
                    && f[3] == row.right_endpoint
                    && f[4] == "spanning_reads",
                "native support template or endpoint differs",
            )?;
            let ids = if f[7] == "." {
                Vec::new()
            } else {
                f[7].split(',').map(str::to_string).collect::<Vec<_>>()
            };
            let count = f[5]
                .parse::<usize>()
                .map_err(|_| err("invalid native support count"))?;
            let seen = per_repeat.entry(row.repeat_node.clone()).or_default();
            require(
                ids.len() == count
                    && ids
                        .iter()
                        .all(|id| !id.is_empty() && seen.insert(id.clone())),
                "native support read IDs are duplicated or disagree with count",
            )?;
            row.path_support = Some(count);
            row.supporting_read_ids = ids;
            row.support_method = f[4].clone();
        }
        Ok(expected)
    }
    fn derived_contract(
        original: &ResolveContract,
        component: &Path,
        checked: &Path,
        details: &Path,
        id: &str,
        selected: &crate::commands::resolve::RepeatPathConstraintResult,
        reverse: bool,
        rotation: usize,
    ) -> Result<ResolveContract, OrgraftError> {
        require(
            selected.search_complete,
            "native repeat candidate search was truncated; source transition remains unresolved",
        )?;
        let text = fs::read_to_string(details)?;
        let mut out = String::new();
        for line in text.lines() {
            let f = line.split('\t').collect::<Vec<_>>();
            let replacement = if f.len() >= 4 && f[1] == id {
                match (f[0], f[2]) {
                    ("resolution", "selected_order") => Some(selected.selected_order.clone()),
                    ("resolution", "oriented_nodes") => Some(selected.oriented_nodes.join(",")),
                    ("resolution", "selected_circular") => Some("true".into()),
                    ("alignment", "reverse_complemented") => Some(reverse.to_string()),
                    ("alignment", "rotation_step") => Some(rotation.to_string()),
                    _ => None,
                }
            } else {
                None
            };
            if let Some(value) = replacement {
                out.push_str(&format!("{}\t{}\t{}\t{}", f[0], f[1], f[2], value));
                for extra in &f[4..] {
                    out.push('\t');
                    out.push_str(extra);
                }
                out.push('\n');
            } else {
                out.push_str(line);
                out.push('\n');
            }
        }
        let contract = contract_from_text(
            &fs::read_to_string(checked)?,
            &fs::read_to_string(component)?,
            &out,
            id,
        )?;
        require(
            contract.path_issues.is_empty() && contract.circular,
            "selected native path is not complete and circular",
        )?;
        require(
            contract.source_sequences == original.source_sequences
                && contract.lineage == original.lineage,
            "source transition changed original nodes or lineage",
        )?;
        let mut walked = contract
            .walk
            .windows(2)
            .map(|w| {
                transition_key(
                    &w[0].source_node,
                    w[0].orientation,
                    &w[1].source_node,
                    w[1].orientation,
                )
            })
            .collect::<Vec<_>>();
        let first = contract.walk.first().unwrap();
        let last = contract.walk.last().unwrap();
        walked.push(transition_key(
            &last.source_node,
            last.orientation,
            &first.source_node,
            first.orientation,
        ));
        walked.sort();
        let mut edges = original
            .topology_links
            .iter()
            .map(physical_key)
            .collect::<Vec<_>>();
        edges.sort();
        require(
            walked == edges,
            "selected ordinary walk omits or duplicates an original physical edge",
        )?;
        Ok(contract)
    }
    // Pure counterpart of the existing node-candidate replay gate. Unlike
    // finish_plan this does not emit diagnostic files during consumer checks.
    fn verify_projection(
        graph: &Gfa,
        target: &FastaRecord,
        contract: &ResolveContract,
    ) -> Result<(), OrgraftError> {
        require(
            graph.segments.keys().collect::<Vec<_>>()
                == contract.source_sequences.keys().collect::<Vec<_>>(),
            "source projection changes node identity",
        )?;
        let sequences = graph
            .segments
            .iter()
            .map(|(n, s)| (n.clone(), s.sequence.to_ascii_uppercase()))
            .collect();
        let edges = contract_edges(graph, contract)?;
        all_edge_overlaps(&sequences, &edges)?;
        let walk = contract
            .walk
            .iter()
            .map(|s| PathStep {
                node: s.source_node.clone(),
                reverse: s.orientation == '-',
            })
            .collect::<Vec<_>>();
        let prepared = prepare_path(&sequences, &edges, &walk, contract.circular).map_err(err)?;
        let target = if contract.alignment_reverse_complemented {
            reverse_complement(&target.sequence)
        } else {
            target.sequence.clone()
        };
        require(
            !target.is_empty() && prepared.source_spell.len() == target.len(),
            "source projection has different full DNA length",
        )?;
        let matches = if contract.circular {
            format!("{target}{target}")
                .find(&prepared.source_spell)
                .is_some_and(|p| p < target.len())
        } else {
            prepared.source_spell == target
        };
        require(
            matches,
            "source projection does not replay the complete original frame",
        )
    }
    fn replay(
        c: &mut Context,
        checked: &Path,
        component: &Path,
        details: &Path,
        id: &str,
        native: &OutputPaths,
    ) -> Result<(ResolveContract, Values), OrgraftError> {
        let mut run = Values::new();
        for f in fields(&native.run_report, 3)? {
            require(
                run.insert(format!("{}/{}", f[0], f[1]), f[2].clone())
                    .is_none(),
                "duplicate rebuild report field",
            )?;
        }
        require(
            get(&run, "run/status")? == "completed"
                && get(&run, "run/subgraph")? == id
                && get(&run, "method/repeat_resolution")?
                    .starts_with(&format!("spanning_read_constrained:{}:", c.repeat))
                && get(&run, "input/checked_gfa")? == "."
                && get(&run, "input/resolve_details")? == ".",
            "source rebuild did not execute the original native resolver",
        )?;
        let first_native = report(&c.first_report)?;
        for (key, path) in [
            ("input/edited_gfa", checked),
            ("input/polished_fasta", c.source.as_path()),
            (
                "input/reads",
                Path::new(get(&first_native, "input/reads/path")?),
            ),
        ] {
            same_path(Path::new(get(&run, key)?), path)?;
        }
        let native_root = native.subgraph_dir.parent().unwrap();
        for (key, path) in [
            ("output/rebuild_gfa", &native.verified_gfa),
            ("output/rebuild_fasta", &native.verified_fasta),
            ("output/repeat_path_support", &native.repeat_path_support),
            ("output/repeat_resolution", &native.repeat_resolution),
        ] {
            let value = PathBuf::from(get(&run, key)?);
            same_path(
                &if value.is_absolute() {
                    value
                } else {
                    native_root.join(value)
                },
                path,
            )?;
        }
        require(
            get(&run, "method/merge_mode")? == "auto_non_repeat_linear_compaction"
                && get(&run, "method/sequence_projection")? != "minimap2_fallback",
            "ordinary transition lacks native verified node projection",
        )?;
        let original = load_contract(checked, component, details, id)?;
        require(
            original.path_issues.is_empty() && original.circular,
            "original path is not legal and circular",
        )?;
        let first = one(&c.source)?;
        let projected = Gfa::read(&native.verified_gfa)?;
        // This is the old input frame. Its VS tags are intentionally not reused.
        verify_projection(&projected, &first, &original)?;
        let source_graph = Gfa::read(component)?;
        let class = original
            .source_sequences
            .keys()
            .map(|n| {
                (
                    n.clone(),
                    if original.walk.iter().filter(|s| &s.source_node == n).count() > 1 {
                        "repeat_node"
                    } else {
                        "single_copy_node"
                    }
                    .into(),
                )
            })
            .collect::<BTreeMap<_, _>>();
        require(
            class.get(&c.repeat).map(String::as_str) == Some("repeat_node"),
            "checkpoint repeat is not an original complete-copy repeat",
        )?;
        // Original reference-aligned draft also replays all original nodes and edges.
        let r = report(&c.first_report)?;
        let draft = one(Path::new(get(&r, "input/draft/path")?))?;
        verify_projection(&source_graph, &draft, &original)?;
        let rows = support(
            &native.repeat_path_support,
            &Gfa::read(component)?,
            &class,
            &native.debug_dir,
        )?;
        let decision = select_dominant_repeat_pairing(&rows, &c.repeat)?;
        let selected = crate::commands::resolve::select_repeat_path_constrained_candidate(
            &crate::commands::resolve::RepeatPathConstraintRequest {
                component_gfa: &native.verified_gfa,
                reference_fasta: &c.source,
                repeat_node: &c.repeat,
                pairing_signature: &decision.pairing_signature,
                max_states: 5000,
                max_candidates: 100,
            },
        )?;
        require(
            selected.search_complete,
            "native repeat candidate search was truncated; source transition remains unresolved",
        )?;
        let (dna, orientation, rotation, _) = normalize_repeat_constrained_sequence(
            &selected.sequence,
            &first.sequence,
            &selected.boundary_offsets,
        )?;
        let normalized = FastaRecord {
            header: first.header.clone(),
            id: first.id.clone(),
            sequence: dna,
        };
        equal_fasta(&normalized, &one(&native.verified_fasta)?)?;
        replay_edits(c, &normalized)?;
        let summary = fields(&native.repeat_resolution, 7)?
            .into_iter()
            .filter(|f| f[0] == "summary")
            .collect::<Vec<_>>();
        require(
            summary.len() == 1
                && summary[0][1] == c.repeat
                && summary[0][2] == decision.pairing_signature
                && summary[0][3] == "Yes",
            "native selection report differs",
        )?;
        let rs = summary[0][6]
            .split(';')
            .filter_map(|s| s.split_once('='))
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect::<Values>();
        for (k, v) in [
            ("method", selected.selection_method.clone()),
            ("local_candidate", selected.selected_candidate.clone()),
            ("selected_order", selected.selected_order.clone()),
            ("oriented_nodes", selected.oriented_nodes.join(",")),
            ("orientation", orientation.to_string()),
            ("rotation_offset", rotation.to_string()),
            ("selected_support", decision.selected_support.to_string()),
            (
                "alternative_support",
                decision.alternative_support.to_string(),
            ),
        ] {
            require(
                get(&rs, k)? == v,
                "native selected operation or frame differs",
            )?;
        }
        for p in [
            &native.verified_gfa,
            &native.verified_fasta,
            &native.repeat_path_support,
            &native.repeat_resolution,
            &native.run_report,
            &native.result_stats,
        ] {
            bind(&mut c.files, p)?;
        }
        bind(&mut c.files, component)?;
        let contract = derived_contract(
            &original,
            component,
            checked,
            details,
            id,
            &selected,
            orientation == '-',
            rotation,
        )?;
        let values = Values::from([
            (
                "schema".into(),
                "orgraft_ordinary_source_transition_v1".into(),
            ),
            ("ordinary_source_dir".into(), c.round.display().to_string()),
            (
                "checked_gfa".into(),
                fs::canonicalize(checked)?.display().to_string(),
            ),
            (
                "component_gfa".into(),
                fs::canonicalize(component)?.display().to_string(),
            ),
            (
                "resolve_details".into(),
                fs::canonicalize(details)?.display().to_string(),
            ),
            ("component".into(), id.into()),
            ("repeat_node".into(), c.repeat.clone()),
            ("selected_order".into(), selected.selected_order),
            ("oriented_nodes".into(), selected.oriented_nodes.join(",")),
            (
                "reverse_complemented".into(),
                (orientation == '-').to_string(),
            ),
            ("rotation_offset".into(), rotation.to_string()),
            (
                "source_frame".into(),
                "original_ordinary_input; not_selected_output_coordinates".into(),
            ),
            (
                "derived_alignment_frame".into(),
                "selected_ordinary_validation_snapshot; final_coordinates_recomputed".into(),
            ),
            (
                "original_alignment_reverse_complemented".into(),
                original.alignment_reverse_complemented.to_string(),
            ),
            (
                "original_alignment_rotation_step".into(),
                original.alignment_rotation_step.to_string(),
            ),
        ]);
        Ok((contract, values))
    }
    fn proof_text(
        values: &Values,
        files: &BTreeMap<PathBuf, String>,
    ) -> Result<String, OrgraftError> {
        let mut text = "metric\tvalue\n".to_string();
        for (k, v) in values {
            require(
                !k.contains(['\t', '\n', '\r']) && !v.contains(['\t', '\n', '\r']),
                "unrepresentable source proof",
            )?;
            text.push_str(&format!("{k}\t{v}\n"));
        }
        for (p, h) in files {
            let p = p.to_str().ok_or_else(|| err("non-UTF8 source"))?;
            require(
                !p.contains(['\t', '\n', '\r']),
                "unrepresentable source path",
            )?;
            text.push_str(&format!("file:{p}\t{h}\n"));
        }
        Ok(text)
    }
    pub(super) fn produce(
        options: &RebuildOptions,
        paths: &OutputPaths,
    ) -> Result<ResolveContract, OrgraftError> {
        let checked = options
            .checked_gfa
            .as_ref()
            .ok_or_else(|| err("checked input is required"))?;
        let details = options
            .resolve_details
            .as_ref()
            .ok_or_else(|| err("resolve input is required"))?;
        let reads = options
            .reads
            .as_ref()
            .ok_or_else(|| err("original reads are required"))?;
        let mut c = context(
            options.ordinary_source_dir.as_ref().unwrap(),
            checked,
            &options.subgraph,
            &options.organelle,
            reads,
            details,
        )?;
        bind(&mut c.files, &options.edited_gfa)?;
        let root = paths.debug_dir.join("ordinary_source");
        let mut nested = options.clone();
        nested.edited_gfa = checked.clone();
        nested.polished_fasta = c.source.clone();
        nested.checked_gfa = None;
        nested.resolve_details = None;
        nested.ordinary_source_dir = None;
        nested.repeat_resolution_node = Some(c.repeat.clone());
        nested.out_dir = root.clone();
        nested.image_reference_fasta = None;
        nested.force = false;
        nested.keep_debug = false;
        stable(&c.files)?;
        run_rebuild(&nested)?;
        let native = OutputPaths::new(&root, &options.subgraph);
        let (contract, mut values) = replay(
            &mut c,
            checked,
            &options.edited_gfa,
            details,
            &options.subgraph,
            &native,
        )?;
        values.insert(
            "reads".into(),
            fs::canonicalize(reads)?.display().to_string(),
        );
        values.insert("organelle".into(), options.organelle.clone());
        stable(&c.files)?;
        let mut out = OpenOptions::new().write(true).create_new(true).open(
            paths
                .logs_dir
                .join(format!("rebuild_{}_ordinary_source.tsv", options.subgraph)),
        )?;
        out.write_all(proof_text(&values, &c.files)?.as_bytes())?;
        out.sync_all()?;
        Ok(contract)
    }
    pub(crate) fn consume(
        proof: &Path,
        checked: &Path,
        component: &Path,
        details: &Path,
        id: &str,
        reads: &Path,
    ) -> Result<ResolveContract, OrgraftError> {
        let before = hash(proof)?;
        let stored = metrics(proof)?;
        require(
            get(&stored, "schema")? == "orgraft_ordinary_source_transition_v1"
                && get(&stored, "component")? == id,
            "wrong ordinary source proof",
        )?;
        for (k, p) in [
            ("checked_gfa", checked),
            ("component_gfa", component),
            ("resolve_details", details),
            ("reads", reads),
        ] {
            same_path(Path::new(get(&stored, k)?), p)?;
        }
        let bound = stored
            .iter()
            .filter_map(|(k, v)| {
                k.strip_prefix("file:")
                    .map(|p| (PathBuf::from(p), v.clone()))
            })
            .collect::<BTreeMap<_, _>>();
        stable(&bound)?;
        let mut c = context(
            Path::new(get(&stored, "ordinary_source_dir")?),
            checked,
            id,
            get(&stored, "organelle")?,
            reads,
            details,
        )?;
        let root = proof
            .parent()
            .and_then(Path::parent)
            .ok_or_else(|| err("proof owner missing"))?;
        let native = OutputPaths::new(&root.join("debug").join(id).join("ordinary_source"), id);
        let (contract, mut values) = replay(&mut c, checked, component, details, id, &native)?;
        values.insert(
            "reads".into(),
            fs::canonicalize(reads)?.display().to_string(),
        );
        values.insert("organelle".into(), get(&stored, "organelle")?.into());
        require(
            fs::read_to_string(proof)? == proof_text(&values, &c.files)?,
            "ordinary source proof does not cover the exact replay inputs",
        )?;
        stable(&c.files)?;
        require(hash(proof)? == before, "source proof changed during replay")?;
        Ok(contract)
    }
    #[cfg(test)]
    mod tests {
        use super::*;
        use std::sync::atomic::{AtomicU64, Ordering};
        static SEQ: AtomicU64 = AtomicU64::new(0);
        struct Fixture {
            root: PathBuf,
            checked: PathBuf,
            component: PathBuf,
            details: PathBuf,
            reads: PathBuf,
            round: PathBuf,
            proof: PathBuf,
            native: OutputPaths,
            checkpoint: PathBuf,
        }
        fn write(path: &Path, text: impl AsRef<[u8]>) {
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, text).unwrap();
        }
        fn dna(mut seed: u64, n: usize) -> String {
            (0..n)
                .map(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    b"ACGT"[((seed >> 32) % 4) as usize] as char
                })
                .collect()
        }
        fn table(rows: &[(&str, String)]) -> String {
            let mut s = "metric\tvalue\n".to_string();
            for (k, v) in rows {
                s.push_str(&format!("{k}\t{v}\n"));
            }
            s
        }
        impl Drop for Fixture {
            fn drop(&mut self) {
                let _ = fs::remove_dir_all(&self.root);
            }
        }
        impl Fixture {
            fn new() -> Self {
                let root = std::env::temp_dir().join(format!(
                    "orgraft-ordinary-source-{}-{}",
                    std::process::id(),
                    SEQ.fetch_add(1, Ordering::Relaxed)
                ));
                fs::create_dir(&root).unwrap();
                let root = fs::canonicalize(root).unwrap();
                let checked =
                    root.join("workflow/mito/subgraph_001/checkpoint_1/checked_draft.gfa");
                let checkpoint = root.join("workflow/mito/subgraph_001/checkpoint_2");
                let component = root.join("resolve/component.gfa");
                let details = root.join("resolve/details.tsv");
                let draft = root.join("resolve/resolved.fasta");
                let reads = root.join("reads.fastq");
                write(&reads, "@fixture\nACGT\n+\nIIII\n");
                let reference = root.join("reference.fasta");
                write(&reference, ">reference\nACGT\n");
                let mut graph = "H\tVN:Z:1.0\n".to_string();
                for (i, n) in ["r0", "a", "b", "r5", "c", "d"].iter().enumerate() {
                    graph.push_str(&format!(
                        "S\t{n}\t{}\n",
                        dna(
                            i as u64 + 7,
                            if n.starts_with('r') {
                                200 + i * 11
                            } else {
                                2200 + i * 17
                            }
                        )
                    ));
                }
                for (a, sa, b, sb) in [
                    ("r0", '+', "c", '+'),
                    ("r0", '+', "d", '+'),
                    ("r0", '-', "a", '+'),
                    ("r0", '-', "b", '+'),
                    ("a", '+', "r5", '-'),
                    ("c", '+', "r5", '+'),
                    ("b", '+', "r5", '+'),
                    ("r5", '+', "d", '-'),
                ] {
                    graph.push_str(&format!("L\t{a}\t{sa}\t{b}\t{sb}\t0M\n"));
                }
                write(&checked, &graph);
                write(&component, &graph);
                let dt=format!("section\tid\tkey\tvalue\tnotes\nrun\t.\tchecked_draft_gfa\t{}\t\nrun\t.\tresolved_fasta\t{}\t\nresolution\tsubgraph_001\tselected_circular\ttrue\t\nresolution\tsubgraph_001\tselected_order\tr5:B->r5_copy1,r0:A->r0_copy1\t\nresolution\tsubgraph_001\toriented_nodes\ta-,r0+,c+,r5+,d-,r0_copy1-,b+,r5_copy1+\t\nresolution\tsubgraph_001\tscore_method\tfixture\t\nalignment\tsubgraph_001\treverse_complemented\tfalse\t\nalignment\tsubgraph_001\trotation_step\t0\t\n",checked.display(),draft.display());
                write(&details, dt);
                let original =
                    load_contract(&checked, &component, &details, "subgraph_001").unwrap();
                let g = Gfa::read(&component).unwrap();
                let walk = original
                    .walk
                    .iter()
                    .map(|s| PathStep {
                        node: s.source_node.clone(),
                        reverse: s.orientation == '-',
                    })
                    .collect::<Vec<_>>();
                let spell = prepare_path(
                    &original.source_sequences,
                    &contract_edges(&g, &original).unwrap(),
                    &walk,
                    true,
                )
                .unwrap()
                .source_spell;
                write(&draft, format!(">subgraph_001\n{spell}\n"));
                let owner = root.join("polish/mito/subgraph_001");
                let source = owner.join("round_1/02.polish/polished_aln.fasta");
                write(&source, format!(">subgraph_001\n{spell}\n"));
                let out = root.join("final");
                let native = OutputPaths::new(
                    &out.join("debug/subgraph_001/ordinary_source"),
                    "subgraph_001",
                );
                native.create(false).unwrap();
                write(&native.verified_gfa, graph);
                let classes = original
                    .source_sequences
                    .keys()
                    .map(|n| {
                        (
                            n.clone(),
                            if n.starts_with('r') {
                                "repeat_node"
                            } else {
                                "single_copy_node"
                            }
                            .into(),
                        )
                    })
                    .collect();
                let mut rows = repeat_path_support_rows(
                    &g,
                    &classes,
                    None,
                    Path::new("unused"),
                    1,
                    &root,
                    false,
                )
                .unwrap();
                for (i, row) in rows.iter_mut().enumerate() {
                    let count = if row.repeat_node == "r5" {
                        if i % 4 == 0 || i % 4 == 3 {
                            4
                        } else {
                            1
                        }
                    } else {
                        2
                    };
                    row.path_support = Some(count);
                    row.path_ratio = Some(0.25);
                    row.support_method = "spanning_reads".into();
                    row.supporting_read_ids = (0..count).map(|n| format!("read_{i}_{n}")).collect();
                }
                write_repeat_path_support(&native.repeat_path_support, &rows).unwrap();
                let decision = select_dominant_repeat_pairing(&rows, "r5").unwrap();
                let selected = crate::commands::resolve::select_repeat_path_constrained_candidate(
                    &crate::commands::resolve::RepeatPathConstraintRequest {
                        component_gfa: &native.verified_gfa,
                        reference_fasta: &source,
                        repeat_node: "r5",
                        pairing_signature: &decision.pairing_signature,
                        max_states: 5000,
                        max_candidates: 100,
                    },
                )
                .unwrap();
                assert!(selected.search_complete);
                let (sequence, orientation, rotation, exact) =
                    normalize_repeat_constrained_sequence(
                        &selected.sequence,
                        &spell,
                        &selected.boundary_offsets,
                    )
                    .unwrap();
                write(
                    &native.verified_fasta,
                    format!(">subgraph_001\n{sequence}\n"),
                );
                let outcome = RepeatResolutionOutcome {
                    decision,
                    selected_candidate: selected.selected_candidate.clone(),
                    local_candidate: selected.selected_candidate.clone(),
                    candidate_count: selected.candidate_count,
                    filtered_candidate_count: selected.filtered_candidate_count,
                    selection_method: selected.selection_method,
                    score: selected.score,
                    continuous_fraction: selected.continuous_fraction,
                    selected_order: selected.selected_order,
                    ordered_nodes: selected.ordered_nodes,
                    oriented_nodes: selected.oriented_nodes,
                    orientation,
                    rotation_offset: rotation,
                    exact_position_matches: exact,
                    candidates: selected.candidates,
                };
                write_repeat_resolution(&native.repeat_resolution, Some(&outcome)).unwrap();
                let nr=format!("section\tkey\tvalue\nrun\tstatus\tcompleted\nrun\tsubgraph\tsubgraph_001\ninput\tchecked_gfa\t.\ninput\tresolve_details\t.\ninput\tedited_gfa\t{}\ninput\tpolished_fasta\t{}\ninput\treads\t{}\nmethod\trepeat_resolution\tspanning_read_constrained:r5:fixture\nmethod\tmerge_mode\tauto_non_repeat_linear_compaction\nmethod\tsequence_projection\tverified_fasta_self_blast_repeat_core\noutput\trebuild_gfa\t{}\noutput\trebuild_fasta\t{}\noutput\trepeat_path_support\t{}\noutput\trepeat_resolution\t{}\n",checked.display(),source.display(),reads.display(),native.verified_gfa.display(),native.verified_fasta.display(),native.repeat_path_support.display(),native.repeat_resolution.display());
                write(&native.run_report, nr);
                write(&native.result_stats, "section\titem\nfixture\tcomplete\n");
                let mismatch =
                    checkpoint.join("round_1/repeat_pairing_repair/repeat_pairing_mismatch.tsv");
                write(
                    &mismatch,
                    "metric\tvalue\nstatus\trepeat_pairing_mismatch\nrepeat_node\tr5\n",
                );
                let old_rebuild = root.join("mutable_old/rebuild.fasta");
                write(
                    &old_rebuild,
                    ">subgraph_001\nTHIS_IS_NOT_THE_ORIGINAL_SNAPSHOT\n",
                );
                let marker = checkpoint.join("repeat_pairing.status.tsv");
                write(
                    &marker,
                    table(&[
                        ("state", "validated".into()),
                        ("repeat_node", "r5".into()),
                        ("source_round", "1".into()),
                        ("validation_round", "3".into()),
                        ("mismatch_report", mismatch.display().to_string()),
                        ("rebuild_fasta", old_rebuild.display().to_string()),
                    ]),
                );
                let corrected_path = checkpoint.join("round_2/polish_aln_v3.fasta");
                let mut corrected = sequence.clone();
                let alternate = if &sequence[9..10] == "A" { "C" } else { "A" };
                corrected.replace_range(9..10, alternate);
                let edits = checkpoint.join("round_2/pos_ref_alt.txt");
                write(&edits, format!("10\t{}\t{alternate}\n", &sequence[9..10]));
                for n in 1..=3 {
                    let dir = owner.join(format!("round_{n}"));
                    let fa = if n == 1 {
                        source.clone()
                    } else {
                        dir.join(format!("01.inputs/linear_subgraph.round_{n}.fasta"))
                    };
                    if n > 1 {
                        write(
                            &fa,
                            format!(
                                ">subgraph_001\n{}\n",
                                if n == 2 { &sequence } else { &corrected }
                            ),
                        );
                    }
                    write(
                        &dir.join("01.inputs/subgraph_reads.fastq.gz"),
                        fs::read(&reads).unwrap(),
                    );
                    write(
                        &dir.join("01.inputs/rotated_reference.fasta"),
                        fs::read(&reference).unwrap(),
                    );
                    let summary = dir.join("03.validate/03.reports/sv_snv_indel_summary.tsv");
                    let high = dir.join("03.validate/03.reports/snv_indel_high.tsv");
                    write(&summary, "metric\tvalue\nstatus\tpass\n");
                    write(
                        &high,
                        if n == 2 {
                            format!("pos\tref\talt\n10\t{}\t{alternate}\n", &sequence[9..10])
                        } else {
                            "pos\tref\talt\n".into()
                        },
                    );
                    let input = if n == 1 {
                        &draft
                    } else if n == 2 {
                        &old_rebuild
                    } else {
                        &corrected_path
                    };
                    let rp=format!("section\tname\tmetric\tvalue\nglobal\trun\torganelle\tmito\nglobal\trun\tsubgraph\tsubgraph_001\nglobal\trun\tvalidate_round\t{n}\nglobal\trun\tper_read_variant_calls\ttrue\noutput\tworkflow_round\tdir\t{}\noutput\tvalidation_fasta\tpath\t{}\ninput\tdraft\tpath\t{}\ninput\treads\tpath\t{}\ninput\treference\tpath\t{}\nstage\tsv-eval:round_{n}\tstatus\tok\nstage\tsnv-indel-eval:round_{n}\tstatus\tok\nsv_eval\tround_1\tsv_support_status\tpass\n",dir.display(),fa.display(),input.display(),reads.display(),reference.display());
                    write(&dir.join("logs/report.tsv"), rp);
                    let mut status = vec![
                        (
                            "status",
                            if n == 1 {
                                "rebuild_ready"
                            } else if n == 2 {
                                "next_round_ready"
                            } else {
                                "complete"
                            }
                            .into(),
                        ),
                        ("snv_indel_high", high.display().to_string()),
                        ("sv_summary", summary.display().to_string()),
                    ];
                    if n == 1 {
                        status.extend([
                            ("sv_correction_mode", "repeat_pairing_rebuild".into()),
                            ("repeat_node", "r5".into()),
                            ("repeat_pairing_source_round", "1".into()),
                            ("repeat_pairing_validation_round", "2".into()),
                            ("repeat_pairing_marker", marker.display().to_string()),
                            (
                                "repeat_pairing_mismatch_report",
                                mismatch.display().to_string(),
                            ),
                            ("corrected_fasta", old_rebuild.display().to_string()),
                        ]);
                    }
                    if n == 2 {
                        status.extend([
                            ("correction_kind", "snv_indel".into()),
                            ("pos_ref_alt", edits.display().to_string()),
                            ("corrected_fasta", corrected_path.display().to_string()),
                        ]);
                    }
                    write(
                        &checkpoint.join(format!("round_{n}/checkpoint_2.status.tsv")),
                        table(&status),
                    );
                }
                let round = owner.join("round_3");
                let proof = out.join("logs/rebuild_subgraph_001_ordinary_source.tsv");
                let f = Self {
                    root,
                    checked,
                    component,
                    details,
                    reads,
                    round,
                    proof,
                    native,
                    checkpoint,
                };
                let mut c = context(
                    &f.round,
                    &f.checked,
                    "subgraph_001",
                    "mito",
                    &f.reads,
                    &f.details,
                )
                .unwrap();
                let (_, mut values) = replay(
                    &mut c,
                    &f.checked,
                    &f.component,
                    &f.details,
                    "subgraph_001",
                    &f.native,
                )
                .unwrap();
                values.insert("reads".into(), f.reads.display().to_string());
                values.insert("organelle".into(), "mito".into());
                write(&f.proof, proof_text(&values, &c.files).unwrap());
                f
            }
            fn consume(&self) -> Result<ResolveContract, OrgraftError> {
                consume(
                    &self.proof,
                    &self.checked,
                    &self.component,
                    &self.details,
                    "subgraph_001",
                    &self.reads,
                )
            }
            fn snapshot(&self) -> BTreeMap<PathBuf, (String, SystemTime)> {
                fn visit(p: &Path, out: &mut BTreeMap<PathBuf, (String, SystemTime)>) {
                    for e in fs::read_dir(p).unwrap() {
                        let p = e.unwrap().path();
                        if p.is_dir() {
                            visit(&p, out);
                        } else {
                            out.insert(
                                p.clone(),
                                (
                                    hash(&p).unwrap(),
                                    fs::metadata(p).unwrap().modified().unwrap(),
                                ),
                            );
                        }
                    }
                }
                let mut out = BTreeMap::new();
                visit(&self.root, &mut out);
                out
            }
            fn resign(&self) {
                let v = metrics(&self.proof).unwrap();
                let mut files = BTreeMap::new();
                let mut values = Values::new();
                for (k, v) in v {
                    if let Some(p) = k.strip_prefix("file:") {
                        let p = PathBuf::from(p);
                        files.insert(p.clone(), hash(&p).unwrap());
                    } else {
                        values.insert(k, v);
                    }
                }
                write(&self.proof, proof_text(&values, &files).unwrap());
            }
        }
        // Optional local integration fixture: original native tables/FASTA, with
        // relocated paths and tiny read/reference stand-ins. This calls no mapper.
        #[test]
        #[ignore = "requires ORGRAFT_ORDINARY_SOURCE_FIXTURE with cached native evidence"]
        fn ordinary_source_replays_cached_native_shape_read_only() {
            let root = PathBuf::from(
                std::env::var("ORGRAFT_ORDINARY_SOURCE_FIXTURE").expect("fixture directory"),
            );
            let source = root.join("original");
            let checked = source.join("workflow/mito/subgraph_001/checkpoint_1/checked_draft.gfa");
            let component =
                source.join("03.resolve_gfa/mito/graph/merged_unresolved_subgraph_001.gfa");
            let details = source.join("03.resolve_gfa/mito/logs/resolve_details.tsv");
            let reads = source.join("01.recruit/mito.fastq.gz");
            let round = source.join("04.polish/mito/subgraph_001/round_3");
            let native = OutputPaths::new(
                &root.join("final/debug/subgraph_001/ordinary_source"),
                "subgraph_001",
            );
            let mut c =
                context(&round, &checked, "subgraph_001", "mito", &reads, &details).unwrap();
            let original = load_contract(&checked, &component, &details, "subgraph_001").unwrap();
            let classes = original
                .source_sequences
                .keys()
                .map(|n| {
                    (
                        n.clone(),
                        if original.walk.iter().filter(|s| &s.source_node == n).count() > 1 {
                            "repeat_node"
                        } else {
                            "single_copy_node"
                        }
                        .into(),
                    )
                })
                .collect();
            let rows = support(
                &native.repeat_path_support,
                &Gfa::read(&component).unwrap(),
                &classes,
                &root,
            )
            .unwrap();
            let decision = select_dominant_repeat_pairing(&rows, &c.repeat).unwrap();
            let selected = crate::commands::resolve::select_repeat_path_constrained_candidate(
                &crate::commands::resolve::RepeatPathConstraintRequest {
                    component_gfa: &native.verified_gfa,
                    reference_fasta: &c.source,
                    repeat_node: &c.repeat,
                    pairing_signature: &decision.pairing_signature,
                    max_states: 5000,
                    max_candidates: 100,
                },
            )
            .unwrap();
            let old = fs::read_to_string(&native.repeat_resolution).unwrap();
            if !old.contains("selected_order=") {
                write(
                    &native.repeat_resolution,
                    old.replace(
                        ";oriented_nodes=",
                        &format!(
                            ";selected_order={};oriented_nodes=",
                            selected.selected_order
                        ),
                    ),
                );
            }
            let (_, mut values) = replay(
                &mut c,
                &checked,
                &component,
                &details,
                "subgraph_001",
                &native,
            )
            .unwrap();
            values.insert("reads".into(), reads.display().to_string());
            values.insert("organelle".into(), "mito".into());
            let proof = root.join("final/logs/rebuild_subgraph_001_ordinary_source.tsv");
            write(&proof, proof_text(&values, &c.files).unwrap());
            fn snapshot(p: &Path, out: &mut BTreeMap<PathBuf, (String, SystemTime)>) {
                for e in fs::read_dir(p).unwrap() {
                    let p = e.unwrap().path();
                    if p.is_dir() {
                        snapshot(&p, out);
                    } else {
                        out.insert(
                            p.clone(),
                            (
                                hash(&p).unwrap(),
                                fs::metadata(p).unwrap().modified().unwrap(),
                            ),
                        );
                    }
                }
            }
            let mut before = BTreeMap::new();
            snapshot(&root, &mut before);
            let derived = consume(
                &proof,
                &checked,
                &component,
                &details,
                "subgraph_001",
                &reads,
            )
            .unwrap();
            assert!(derived.path_issues.is_empty());
            assert_eq!(derived.walk.len(), 8);
            let mut after = BTreeMap::new();
            snapshot(&root, &mut after);
            assert_eq!(before, after);
        }
        #[test]
        fn ordinary_source_does_not_consult_stale_history_when_original_path_passes() {
            let f = Fixture::new();
            let source = f
                .round
                .parent()
                .unwrap()
                .join("round_1/02.polish/polished_aln.fasta");
            let args = vec![
                "--edited-gfa".into(),
                f.component.display().to_string(),
                "--checked-gfa".into(),
                f.checked.display().to_string(),
                "--resolve-details".into(),
                f.details.display().to_string(),
                "--polished-fasta".into(),
                source.display().to_string(),
                "--reads".into(),
                f.reads.display().to_string(),
                "--ordinary-source-dir".into(),
                f.root.join("invalid/stale/round_9").display().to_string(),
            ];
            let options = RebuildOptions::from_args(&args).unwrap();
            let mut contract =
                load_contract(&f.checked, &f.component, &f.details, "subgraph_001").unwrap();
            let mut graph = Gfa::read(&f.component).unwrap();
            let target = one(&source).unwrap();
            let paths = OutputPaths::new(&f.root.join("baseline-success"), "subgraph_001");
            let result = recorded_candidate_or_source(
                &mut graph,
                &target,
                &source,
                &mut contract,
                Path::new("/not/a/mapper"),
                Path::new("/not/blast"),
                &options,
                &paths,
            )
            .unwrap();
            assert_eq!(result.0.sequences, contract.source_sequences);
            assert!(!paths.debug_dir.join("ordinary_source").exists());
            assert!(!paths
                .logs_dir
                .join("rebuild_subgraph_001_ordinary_source.tsv")
                .exists());
        }
        #[test]
        fn ordinary_source_replays_selected_transition_and_edits_read_only() {
            let f = Fixture::new();
            let before = f.snapshot();
            let contract = f.consume().unwrap();
            assert!(contract.path_issues.is_empty());
            assert_eq!(before, f.snapshot());
            assert!(!f.native.debug_dir.exists());
            assert!(fs::read_to_string(&f.proof)
                .unwrap()
                .contains("not_selected_output_coordinates"));
            assert!(!fs::read_to_string(&f.proof)
                .unwrap()
                .contains("mutable_old/rebuild.fasta\t"));
        }
        #[test]
        fn ordinary_source_rejects_changed_graph_marker_reads_and_edit_chain_without_writes() {
            for name in [
                "checked",
                "marker",
                "reads",
                "selected",
                "ref",
                "frame",
                "missing_edge",
            ] {
                let f = Fixture::new();
                let p = match name {
                    "checked" => f.checked.clone(),
                    "marker" => f.checkpoint.join("repeat_pairing.status.tsv"),
                    "reads" => f.reads.clone(),
                    "selected" => f.native.repeat_resolution.clone(),
                    "ref" => f.checkpoint.join("round_2/pos_ref_alt.txt"),
                    "frame" => f.round.join("01.inputs/linear_subgraph.round_3.fasta"),
                    _ => f.native.verified_gfa.clone(),
                };
                let text = fs::read_to_string(&p).unwrap();
                let changed = match name {
                    "checked" => text.replacen("\t0M", "\t1M", 1),
                    "marker" => text.replace("repeat_node\tr5", "repeat_node\tr0"),
                    "reads" => text.replace("ACGT", "TGCA"),
                    "selected" => text.replace("selected_order=r5:", "selected_order=r0:"),
                    "ref" => text.replacen("10\t", "11\t", 1),
                    "frame" => text.replace("\n", "\nA"),
                    _ => text
                        .lines()
                        .filter(|l| !l.starts_with("L\tr0\t+\tc"))
                        .collect::<Vec<_>>()
                        .join("\n"),
                };
                assert_ne!(text, changed, "{name}");
                write(&p, changed);
                f.resign();
                let before = f.snapshot();
                assert!(f.consume().is_err(), "{name}");
                assert_eq!(before, f.snapshot(), "consumer wrote on {name} failure");
            }
        }
        #[test]
        fn ordinary_source_requires_complete_candidate_search_and_exact_full_frame() {
            let f = Fixture::new();
            let result = crate::commands::resolve::select_repeat_path_constrained_candidate(
                &crate::commands::resolve::RepeatPathConstraintRequest {
                    component_gfa: &f.native.verified_gfa,
                    reference_fasta: &f
                        .round
                        .parent()
                        .unwrap()
                        .join("round_1/02.polish/polished_aln.fasta"),
                    repeat_node: "r5",
                    pairing_signature: "unused",
                    max_states: 1,
                    max_candidates: 1,
                },
            );
            assert!(result.is_err());
            let raw = Gfa::read(&f.component).unwrap();
            let classes = raw
                .segments
                .keys()
                .map(|n| {
                    (
                        n.clone(),
                        if n.starts_with('r') {
                            "repeat_node"
                        } else {
                            "single_copy_node"
                        }
                        .into(),
                    )
                })
                .collect();
            let rows = support(&f.native.repeat_path_support, &raw, &classes, &f.root).unwrap();
            let target = rows
                .iter()
                .filter(|r| r.repeat_node == "r5")
                .collect::<Vec<_>>();
            let mut found = false;
            for (a, b) in [(0, 3), (1, 2)] {
                let mut signatures = [
                    repeat_path_physical_signature(target[a]),
                    repeat_path_physical_signature(target[b]),
                ];
                signatures.sort();
                for states in 2..30 {
                    if let Ok(selected) =
                        crate::commands::resolve::select_repeat_path_constrained_candidate(
                            &crate::commands::resolve::RepeatPathConstraintRequest {
                                component_gfa: &f.native.verified_gfa,
                                reference_fasta: &f
                                    .round
                                    .parent()
                                    .unwrap()
                                    .join("round_1/02.polish/polished_aln.fasta"),
                                repeat_node: "r5",
                                pairing_signature: &signatures.join("|"),
                                max_states: states,
                                max_candidates: 1,
                            },
                        )
                    {
                        if !selected.search_complete {
                            let base =
                                load_contract(&f.checked, &f.component, &f.details, "subgraph_001")
                                    .unwrap();
                            assert!(derived_contract(
                                &base,
                                &f.component,
                                &f.checked,
                                &f.details,
                                "subgraph_001",
                                &selected,
                                false,
                                0
                            )
                            .is_err());
                            found = true;
                            break;
                        }
                    }
                }
            }
            assert!(
                found,
                "fixture must exercise a nonempty truncated candidate set"
            );
            let original =
                load_contract(&f.checked, &f.component, &f.details, "subgraph_001").unwrap();
            let g = Gfa::read(&f.component).unwrap();
            let mut target = one(&f
                .round
                .parent()
                .unwrap()
                .join("round_1/02.polish/polished_aln.fasta"))
            .unwrap();
            verify_projection(&g, &target, &original).unwrap();
            target.sequence.push('A');
            assert!(verify_projection(&g, &target, &original).is_err());
        }
    }
}

pub(crate) fn ordinary_source_checkpoint(checked: &Path) -> Result<Option<PathBuf>,OrgraftError> {
    let checked=fs::canonicalize(checked)?;
    if checked.file_name().and_then(|n|n.to_str())!=Some("checked_draft.gfa") || checked.parent().and_then(Path::file_name).and_then(|n|n.to_str())!=Some("checkpoint_1") { return Ok(None); }
    Ok(Some(ordinary_source::checkpoint(&checked)?))
}
pub(crate) fn load_ordinary_source_contract(proof: &Path, checked: &Path, component: &Path, details: &Path, id: &str, reads: &Path) -> Result<ResolveContract,OrgraftError> {
    ordinary_source::consume(proof,checked,component,details,id,reads)
}

#[cfg(test)]
mod tests_multi_context {
    use super::*;
    fn templates() -> Vec<RepeatPathTemplate> {
        (0..3).map(|i|RepeatPathTemplate{id:format!("c{i}"),repeat_node:"R".into(),row_index:if i<2 {0}else{1},
            sequence:"A".repeat(4189),required_start:1000,required_end:3189}).collect()
    }
    fn paf(read:&str,target:&str,matches:usize,start:usize,end:usize,primary:bool)->String {
        format!("{read}\t5000\t0\t{}\t+\t{target}\t4189\t{start}\t{end}\t{matches}\t{}\t60\ttp:A:{}\n",end-start,end-start,if primary{'P'}else{'S'})
    }
    #[test]
    fn multi_context_same_physical_row_is_counted_once_with_original_rank() {
        let t=templates();let p=paf("read","c0",4000,0,4000,true)+&paf("read","c1",4000,0,4000,false);
        let got=repeat_path_support_from_paf(&t,&p,true).unwrap();assert_eq!(got.values().map(Vec::len).sum::<usize>(),1);
        assert_eq!(got["c0"],vec!["read"]);
        assert!(repeat_path_support_from_paf(&t,&p,false).unwrap().is_empty());
        let p=p.clone()+&p;assert_eq!(repeat_path_support_from_paf(&t,&p,true).unwrap(),got);
    }
    #[test]
    fn multi_context_other_physical_row_keeps_original_ambiguity_and_windows() {
        let t=templates();let p=paf("read","c0",4000,0,4000,true)+&paf("read","c1",4000,0,4000,false)+&paf("read","c2",4000,0,4000,true);
        assert!(repeat_path_support_from_paf(&t,&p,true).unwrap().is_empty());
        let p=paf("low_identity","c0",2000,0,4000,true)+&paf("left_short","c0",2800,1200,4000,true)+&paf("right_short","c0",3000,0,3000,true);
        assert!(repeat_path_support_from_paf(&t,&p,true).unwrap().is_empty());
        let p=paf("one","c0",4000,0,4000,true)+&paf("two","c2",4000,0,4000,true);
        assert_eq!(repeat_path_support_from_paf(&t,&p,true).unwrap(),repeat_path_support_from_paf(&t,&p,false).unwrap());
    }
    #[test]
    fn multi_context_secondary_saturation_rejects_the_whole_measurement() {
        let t=templates();
        let primary=paf("read","c0",4000,0,4000,true);
        let secondary=paf("read","c1",4000,0,4000,false);
        let below=primary.clone()+&secondary.repeat(REPEAT_PATH_SECONDARY_LIMIT-1);
        assert_eq!(repeat_path_support_from_paf(&t,&below,true).unwrap()["c0"],vec!["read"]);
        let saturated=below+&secondary;
        assert!(repeat_path_support_from_paf(&t,&saturated,true).unwrap_err().to_string().contains("saturated secondary"));
        assert!(repeat_path_support_from_paf(&t,&saturated,false).is_ok());
        // Check the raw mapper output before identity, window, or known-query
        // filters: those discarded alignments can still consume the -N cap.
        let censored=primary+&paf("read","unknown",1,1200,1202,false).repeat(REPEAT_PATH_SECONDARY_LIMIT);
        assert!(repeat_path_support_from_paf(&t,&censored,true).is_err());
        let many_reads=(0..REPEAT_PATH_SECONDARY_LIMIT+1).map(|i|paf(&format!("r{i}"),"c0",4000,0,4000,false)).collect::<String>();
        assert_eq!(repeat_path_support_from_paf(&t,&many_reads,true).unwrap()["c0"].len(),REPEAT_PATH_SECONDARY_LIMIT+1);
    }

}
