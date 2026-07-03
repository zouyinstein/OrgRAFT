# source /Users/zouyinstein-home/Documents/Codex/GFA_Editor/.venv/bin/activate
# ln -sf /Users/zouyinstein-home/Documents/Codex/OrgRAFT soft_paths.txt
# export PATH=/Users/zouyinstein-home/Documents/Codex/OrgRAFT:$PATH
# orgraft setup --soft-paths soft_paths.txt --requirements  /Users/zouyinstein-home/Documents/Codex/OrgRAFT/requirements.txt

# mkdir results_workflow
# orgraft workflow init --out results_workflow/orgraft.workflow.toml
# 修改reads,baits
# orgraft workflow plan --config results_workflow/orgraft.workflow.toml
# bash results_workflow/workflow.commands.sh
# 修改文件/Users/zouyinstein-home/Documents/Codex/OrgRAFT/results_workflow/02.draft_asm/mito/03.finalize_graph/graph.edited.gfa
# 修改results_workflow/orgraft.workflow.toml
#   draft_graph = "${results_dir}/draft_asm/${organelle}/03.finalize_graph/graph.edited.gfa"


# orgraft workflow checkpoint1 --config results_workflow/orgraft.workflow.toml --case mito_subgraph_001 --force
# orgraft workflow plan --config results_workflow/orgraft.workflow.toml --case mito_subgraph_001
# bash /Users/zouyinstein-home/Documents/Codex/OrgRAFT/results_workflow/workflow/mito/subgraph_001/workflow.commands.sh
# orgraft workflow runtime-summary --config results_workflow/orgraft.workflow.toml --force
