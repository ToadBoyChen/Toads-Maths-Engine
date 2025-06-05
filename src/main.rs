use docker_rust_bin::*;

fn main() {
    let group = IntegerAddition;
    let a = 5;
    let b = 3;
    println!("{} + {} = {}", a, b, group.operate_binary(a, b));
    println!("Identity: {}", group.identity());
    println!("Inverse of {}: {}", a, group.inverse(a));
}
