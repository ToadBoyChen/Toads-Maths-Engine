use docker_rust_bin::*;

#[test]
fn test_mod7_ideal() {
    // Test the trivial ideal {0}
    let ideal = [0];
    let mod7 = Mod7;
    for &a in &ideal {
        for &b in &ideal {
            assert!(
                ideal.contains(&mod7.operate_binary(a, b)),
                "Ideal not closed under addition: {} + {} mod 7 = {}",
                a, b, mod7.operate_binary(a, b)
            );
        }
    }
    for r in 0..7 {
        for &i in &ideal {
            let prod = (r * i) % 7;
            assert!(
                ideal.contains(&prod),
                "Ideal not closed under multiplication: {} * {} mod 7 = {}",
                r, i, prod
            );
        }
    }
}

#[test]
fn test_mod7_non_ideal() {
    let not_ideal = [0, 2, 4];
    let mod7 = Mod7;
    let mut fails = false;
    'outer: for r in 0..7 {
        for &i in &not_ideal {
            let prod = (r * i) % 7;
            if !not_ideal.contains(&prod) {
                fails = true;
                break 'outer;
            }
        }
    }
    assert!(fails, "{{0,2,4}} should not be an ideal in Z/7Z");
}
#[test]
fn test_mod7_ideal_edge_cases() {
    let mod7 = Mod7;
    let empty: [i32; 0] = [];
    for &a in &empty {
        for &b in &empty {
            assert!(empty.contains(&mod7.operate_binary(a, b)));
        }
    }
    let full: Vec<i32> = (0..7).collect();
    for &a in &full {
        for &b in &full {
            assert!(full.contains(&mod7.operate_binary(a, b)));
        }
    }
}
