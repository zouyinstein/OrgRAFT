pub struct CommandContract {
    pub command: &'static str,
    pub origin: &'static str,
    pub purpose: &'static str,
    pub inputs: &'static [&'static str],
    pub outputs: &'static [&'static str],
    pub notes: &'static [&'static str],
}

pub fn print_contract(contract: &CommandContract) {
    println!("orgraft {}", contract.command);
    println!();
    println!("Origin: {}", contract.origin);
    println!("Purpose: {}", contract.purpose);
    println!();
    print_list("Inputs", contract.inputs);
    print_list("Outputs", contract.outputs);
    print_list("Notes", contract.notes);
}

fn print_list(label: &str, values: &[&str]) {
    if values.is_empty() {
        return;
    }

    println!("{label}:");
    for value in values {
        println!("  - {value}");
    }
    println!();
}
