#[test]
fn test_modn_ideal() {
    let ideal = [0]; // the trivial ideal {0}
    let modn = ModN::new(7);

    for &a in &ideal {
        for &b in &ideal {
            let sum = modn.operate_binary(a, b);
            assert!(
                ideal.contains(&sum),
                "Ideal not closed under addition: {} + {} mod {} = {}",
                a, b, modn.n, sum
            );
        }
    }

    for r in 0..modn.n {
        for &i in &ideal {
            let prod = (r * i) % modn.n;
            assert!(
                ideal.contains(&prod),
                "Ideal not closed under multiplication: {} * {} mod {} = {}",
                r, i, modn.n, prod
            );
        }
    }
}

#[test]
fn test_modn_non_ideal() {
    let not_ideal = [0, 2, 4]; // not an ideal in Z/7Z
    let modn = ModN::new(7);
    let mut fails = false;

    'outer: for r in 0..modn.n {
        for &i in &not_ideal {
            let prod = (r * i) % modn.n;
            if !not_ideal.contains(&prod) {
                fails = true;
                break 'outer;
            }
        }
    }

    assert!(
        fails,
        "{{0,2,4}} should not be an ideal in Z/{}Z", modn.n
    );
}

#[test]
fn test_modn_ideal_edge_cases() {
    let modn = ModN::new(7);
    
    let empty: [i32; 0] = [];
    for &a in &empty {
        for &b in &empty {
            assert!(empty.contains(&modn.operate_binary(a, b)));
        }
    }

    let full: Vec<i32> = (0..modn.n).collect();
    for &a in &full {
        for &b in &full {
            assert!(
                full.contains(&modn.operate_binary(a, b)),
                "Full set should be closed under addition mod {}: {} + {} = {}",
                modn.n, a, b, modn.operate_binary(a, b)
            );
        }
    }
}
