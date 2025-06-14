#[cfg(test)]
mod tests {
    include!("tests/linear_algebra_tests.rs");
    include!("tests/number_system_tests.rs");
    include!("tests/validating_structures_tests.rs");
    include!("tests/p2p_tests.rs");
}

mod libs {
    include!("lib/algebras.rs");
    include!("lib/specialised_algebras.rs");
}

use crate::libs::*;

fn main() {
    let group = IntegerAddition;
    let a = 5;
    let b = 3;
    println!("{} + {} = {}", a, b, group.operate_binary(a, b));
    println!("Identity: {}", group.identity());
    println!("Inverse of {}: {}", a, group.inverse(a));
}