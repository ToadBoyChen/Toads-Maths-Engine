use docker_rust_bin::*;
use rand::Rng;

// IntegerAddition implements Group<i32>
#[test]
fn test_integer_addition_group_laws() {
    let add = IntegerAddition;

    // Closure
    for a in -100..100 {
        for b in -100..100 {
            let res = add.operate_binary(a, b);
            // In Z, always closed
            assert!(res >= i32::MIN && res <= i32::MAX);
        }
    }

    // Associativity
    for a in -10..10 {
        for b in -10..10 {
            for c in -10..10 {
                assert_eq!(
                    add.operate_binary(add.operate_binary(a, b), c),
                    add.operate_binary(a, add.operate_binary(b, c))
                );
            }
        }
    }

    // Commutativity
    for a in -100..100 {
        for b in -100..100 {
            assert_eq!(add.operate_binary(a, b), add.operate_binary(b, a));
        }
    }

    // Identity
    assert_eq!(add.identity(), 0);
    for a in -100..100 {
        assert_eq!(add.operate_binary(a, add.identity()), a);
        assert_eq!(add.operate_binary(add.identity(), a), a);
    }

    // Inverse
    for a in -100..100 {
        assert_eq!(add.operate_binary(a, add.inverse(a)), add.identity());
    }

    // Edge cases: empty and single-element slices
    assert_eq!(add.operate(&[]), 0);
    assert_eq!(add.operate(&[42]), 42);

    // Large numbers
    let big = i32::MAX - 1;
    assert_eq!(add.operate_binary(big, 1), i32::MAX);
}

// IntegerMultiplication is a monoid, not a group (no inverse for 0)
#[test]
fn test_integer_multiplication_monoid_laws() {
    let mul = IntegerMultiplication;

    // Closure
    for a in -100..100 {
        for b in -100..100 {
            let res = mul.operate_binary(a, b);
            assert!(res >= i32::MIN && res <= i32::MAX);
        }
    }

    // Associativity
    for a in -10..10 {
        for b in -10..10 {
            for c in -10..10 {
                assert_eq!(
                    mul.operate_binary(mul.operate_binary(a, b), c),
                    mul.operate_binary(a, mul.operate_binary(b, c))
                );
            }
        }
    }

    // Commutativity
    for a in -100..100 {
        for b in -100..100 {
            assert_eq!(mul.operate_binary(a, b), mul.operate_binary(b, a));
        }
    }

    // Identity
    assert_eq!(mul.operate_binary(1, 1), 1);
    for a in -100..100 {
        assert_eq!(mul.operate_binary(a, 1), a);
        assert_eq!(mul.operate_binary(1, a), a);
    }

    // Edge cases
    assert_eq!(mul.operate(&[]), 1);
    assert_eq!(mul.operate(&[42]), 42);

    // Large numbers
    let big = i32::MAX / 2;
    assert_eq!(mul.operate_binary(big, 2), big * 2);
}

// Modular arithmetic (mod 7) as a group under addition
#[test]
fn test_mod7_group_laws() {
    let mod7 = Mod7;

    // Closure
    for a in 0..7 {
        for b in 0..7 {
            let res = mod7.operate_binary(a, b);
            assert!((0..7).contains(&res));
        }
    }

    // Associativity
    for a in 0..7 {
        for b in 0..7 {
            for c in 0..7 {
                assert_eq!(
                    mod7.operate_binary(mod7.operate_binary(a, b), c),
                    mod7.operate_binary(a, mod7.operate_binary(b, c))
                );
            }
        }
    }

    // Commutativity
    for a in 0..7 {
        for b in 0..7 {
            assert_eq!(mod7.operate_binary(a, b), mod7.operate_binary(b, a));
        }
    }

    // Identity
    for a in 0..7 {
        assert_eq!(mod7.operate_binary(a, 0), a);
        assert_eq!(mod7.operate_binary(0, a), a);
    }

    // Inverse
    for a in 0..7 {
        let inv = (7 - a) % 7;
        assert_eq!(mod7.operate_binary(a, inv), 0);
    }

    // Edge cases
    assert_eq!(mod7.operate(&[]), 0);
    assert_eq!(mod7.operate(&[3]), 3);
}

// Complex addition: group laws (closure, associativity, commutativity, identity, inverse)
#[test]
fn test_complex_addition_group_laws() {
    let add = ComplexAddition;
    let mut rng = rand::thread_rng();

    // Randomized closure, associativity, commutativity, identity, inverse
    for _ in 0..100 {
        let a = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        let b = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        let c = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };

        // Closure
        let res = add.operate_binary(a, b);
        // Always a Complex

        // Associativity
        assert_eq!(
            add.operate_binary(add.operate_binary(a, b), c),
            add.operate_binary(a, add.operate_binary(b, c))
        );

        // Commutativity
        assert_eq!(add.operate_binary(a, b), add.operate_binary(b, a));

        // Identity
        let zero = Complex { re: 0.0, im: 0.0 };
        assert_eq!(add.operate_binary(a, zero), a);
        assert_eq!(add.operate_binary(zero, a), a);

        // Inverse
        let inv = Complex { re: -a.re, im: -a.im };
        assert_eq!(add.operate_binary(a, inv), zero);
    }

    // Edge cases
    let zero = Complex { re: 0.0, im: 0.0 };
    assert_eq!(add.operate(&[]), zero);
    let a = Complex { re: 1.0, im: 2.0 };
    assert_eq!(add.operate(&[a]), a);
}

// Matrix addition: group laws (closure, associativity, commutativity, identity, inverse)
#[test]
fn test_matrix_addition_group_laws() {
    let add = MatrixAddition;
    let a = Matrix2x2([[1, 2], [3, 4]]);
    let b = Matrix2x2([[5, 6], [7, 8]]);
    let c = Matrix2x2([[9, 10], [11, 12]]);
    let zero = Matrix2x2([[0, 0], [0, 0]]);

    // Closure
    let res = add.operate_binary(a, b);
    // Always a Matrix2x2

    // Associativity
    assert_eq!(
        add.operate_binary(add.operate_binary(a, b), c),
        add.operate_binary(a, add.operate_binary(b, c))
    );

    // Commutativity
    assert_eq!(add.operate_binary(a, b), add.operate_binary(b, a));

    // Identity
    assert_eq!(add.operate_binary(a, zero), a);
    assert_eq!(add.operate_binary(zero, a), a);

    // Inverse
    let inv = Matrix2x2([[-1, -2], [-3, -4]]);
    assert_eq!(add.operate_binary(a, inv), zero);

    // Edge cases
    assert_eq!(add.operate(&[]), zero);
    assert_eq!(add.operate(&[a]), a);
}
