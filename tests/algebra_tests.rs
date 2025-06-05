use docker_rust_bin::*;
use rand::Rng;

macro_rules! assert_complex_eq {
    ($a:expr, $b:expr) => {
        assert!(($a.re - $b.re).abs() < 1e-6 && ($a.im - $b.im).abs() < 1e-6,
            "Complex mismatch: {:?} != {:?}", $a, $b);
    };
}

// IntegerAddition implements Group<i32>
#[test]
fn test_integer_addition_group_laws() {
    let add = IntegerAddition;

    // Closure + edge case protection
    for &a in &[i32::MIN, -100, -1, 0, 1, 100, i32::MAX] {
        for &b in &[i32::MIN, -100, -1, 0, 1, 100, i32::MAX] {
            let res = add.operate_binary(a, b);
            // Rust wraps on overflow, but in abstract algebra, this should be closed in Z
            assert_eq!(res, a.wrapping_add(b));
        }
    }

    // Associativity
    for a in -20..20 {
        for b in -20..20 {
            for c in -20..20 {
                let ab = add.operate_binary(a, b);
                let bc = add.operate_binary(b, c);
                assert_eq!(
                    add.operate_binary(ab, c),
                    add.operate_binary(a, bc)
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
    for a in -200..200 {
        assert_eq!(add.operate_binary(a, 0), a);
        assert_eq!(add.operate_binary(0, a), a);
    }

    // Inverse
    for a in -100..100 {
        assert_eq!(add.operate_binary(a, -a), 0);
    }

    // Operate slice
    assert_eq!(add.operate(&[]), 0);
    assert_eq!(add.operate(&[42]), 42);

    // Large value sanity
    assert_eq!(add.operate_binary(i32::MAX - 1, 1), i32::MAX);
}

// IntegerMultiplication is a commutative monoid, not a group
#[test]
fn test_integer_multiplication_monoid_laws() {
    let mul = IntegerMultiplication;

    // Closure
    for &a in &[i32::MIN, -100, -1, 0, 1, 100, i32::MAX] {
        for &b in &[i32::MIN, -100, -1, 0, 1, 100, i32::MAX] {
            let res = mul.operate_binary(a, b);
            assert_eq!(res, a.wrapping_mul(b));
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
    for a in -200..200 {
        assert_eq!(mul.operate_binary(a, 1), a);
        assert_eq!(mul.operate_binary(1, a), a);
    }

    // Edge
    assert_eq!(mul.operate(&[]), 1);
    assert_eq!(mul.operate(&[7]), 7);

    let big = i32::MAX / 2;
    assert_eq!(mul.operate_binary(big, 2), big * 2);
}

// Complex addition: use approx_eq
#[test]
fn test_complex_addition_group_laws() {
    let add = ComplexAddition;
    let mut rng = rand::thread_rng();
    let zero = Complex { re: 0.0, im: 0.0 };

    for _ in 0..100 {
        let a = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        let b = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        let c = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };

        // Associativity
        assert_complex_eq!(
            add.operate_binary(add.operate_binary(a, b), c),
            add.operate_binary(a, add.operate_binary(b, c))
        );

        // Commutativity
        assert_complex_eq!(add.operate_binary(a, b), add.operate_binary(b, a));

        // Identity
        assert_complex_eq!(add.operate_binary(a, zero), a);
        assert_complex_eq!(add.operate_binary(zero, a), a);

        // Inverse
        let inv = Complex { re: -a.re, im: -a.im };
        assert_complex_eq!(add.operate_binary(a, inv), zero);
    }

    assert_complex_eq!(add.operate(&[]), zero);
    let a = Complex { re: 5.0, im: 7.0 };
    assert_complex_eq!(add.operate(&[a]), a);
}

// Matrix group: integer matrix addition
#[test]
fn test_matrix_addition_group_laws() {
    let add = MatrixAddition;
    let a = Matrix2x2([[1, 2], [3, 4]]);
    let b = Matrix2x2([[5, 6], [7, 8]]);
    let c = Matrix2x2([[9, 10], [11, 12]]);
    let zero = Matrix2x2([[0, 0], [0, 0]]);

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

    // Operate slice
    assert_eq!(add.operate(&[]), zero);
    assert_eq!(add.operate(&[a]), a);
}

#[test]
fn test_modn_group_laws_for_n7() {
    let modn = ModN::new(7);

    for a in 0..modn.n {
        for b in 0..modn.n {
            let res = modn.operate_binary(a, b);
            assert!((0..modn.n).contains(&res));

            // Commutativity
            assert_eq!(modn.operate_binary(a, b), modn.operate_binary(b, a));
        }
    }

    for a in 0..modn.n {
        for b in 0..modn.n {
            for c in 0..modn.n {
                let ab = modn.operate_binary(a, b);
                let bc = modn.operate_binary(b, c);
                assert_eq!(
                    modn.operate_binary(ab, c),
                    modn.operate_binary(a, bc)
                );
            }
        }

        // Identity
        assert_eq!(modn.operate_binary(a, 0), a);
        assert_eq!(modn.operate_binary(0, a), a);

        // Inverse
        let inv = modn.inverse(a);
        assert_eq!(modn.operate_binary(a, inv), 0);
    }

    assert_eq!(modn.operate(&[]), 0);
    assert_eq!(modn.operate(&[6]), 6);
}

#[test]
fn test_modn_small_groups() {
    for n in 1..=3 {
        let modn = ModN::new(n);

        for a in 0..n {
            for b in 0..n {
                // Closure
                let ab = modn.operate_binary(a, b);
                assert!((0..n).contains(&ab));

                // Commutativity
                assert_eq!(ab, modn.operate_binary(b, a));

                // Associativity
                let ab_c = modn.operate_binary(ab, 0);
                let a_bc = modn.operate_binary(a, modn.operate_binary(b, 0));
                assert_eq!(ab_c, a_bc);
            }

            // Identity
            assert_eq!(modn.operate_binary(a, modn.identity()), a);

            // Inverse
            let inv = modn.inverse(a);
            assert_eq!(modn.operate_binary(a, inv), 0);
        }
    }
}

#[test]
#[should_panic(expected = "Modulus must be positive and non-zero")]
fn test_mod0_should_panic() {
    let _ = ModN::new(0);
}

#[test]
fn test_quaternion_addition_group_laws() {
    let add = QuaternionAddition;
    let zero = Quaternion { re: 0.0, i: 0.0, j: 0.0, k: 0.0 };
    let mut rng = rand::thread_rng();

    for _ in 0..100 {
        let a = Quaternion {
            re: rng.gen_range(-100.0..100.0),
            i: rng.gen_range(-100.0..100.0),
            j: rng.gen_range(-100.0..100.0),
            k: rng.gen_range(-100.0..100.0),
        };
        let b = Quaternion {
            re: rng.gen_range(-100.0..100.0),
            i: rng.gen_range(-100.0..100.0),
            j: rng.gen_range(-100.0..100.0),
            k: rng.gen_range(-100.0..100.0),
        };
        let c = Quaternion {
            re: rng.gen_range(-100.0..100.0),
            i: rng.gen_range(-100.0..100.0),
            j: rng.gen_range(-100.0..100.0),
            k: rng.gen_range(-100.0..100.0),
        };

        // Associativity
        let ab = add.operate_binary(a, b);
        let bc = add.operate_binary(b, c);
        let ab_c = add.operate_binary(ab, c);
        let a_bc = add.operate_binary(a, bc);
        assert_eq!(ab_c, a_bc);

        // Commutativity
        assert_eq!(add.operate_binary(a, b), add.operate_binary(b, a));

        // Identity
        assert_eq!(add.operate_binary(a, zero), a);
        assert_eq!(add.operate_binary(zero, a), a);

        // Inverse
        let inv = Quaternion { re: -a.re, i: -a.i, j: -a.j, k: -a.k };
        assert_eq!(add.operate_binary(a, inv), zero);
    }

    // Closure
    let a = Quaternion { re: 1.0, i: 2.0, j: 3.0, k: 4.0 };
    let b = Quaternion { re: -1.0, i: -2.0, j: -3.0, k: -4.0 };
    let result = add.operate_binary(a, b);
    assert_eq!(result, zero);

    // Slice ops
    assert_eq!(add.operate(&[]), zero);
    assert_eq!(add.operate(&[a]), a);
    assert_eq!(add.operate(&[a, b]), zero);
}

#[test]
fn test_quaternion_multiplication_known_identities() {
    let mul = QuaternionMultiplication;

    let one = Quaternion { re: 1.0, i: 0.0, j: 0.0, k: 0.0 };
    let i = Quaternion { re: 0.0, i: 1.0, j: 0.0, k: 0.0 };
    let j = Quaternion { re: 0.0, i: 0.0, j: 1.0, k: 0.0 };
    let k = Quaternion { re: 0.0, i: 0.0, j: 0.0, k: 1.0 };
    let neg_one = Quaternion { re: -1.0, i: 0.0, j: 0.0, k: 0.0 };

    // Identity
    assert_eq!(mul.operate_binary(one, i), i);
    assert_eq!(mul.operate_binary(i, one), i);

    assert_eq!(mul.operate_binary(one, j), j);
    assert_eq!(mul.operate_binary(j, one), j);

    assert_eq!(mul.operate_binary(one, k), k);
    assert_eq!(mul.operate_binary(k, one), k);

    // Squares
    assert_eq!(mul.operate_binary(i, i), neg_one);
    assert_eq!(mul.operate_binary(j, j), neg_one);
    assert_eq!(mul.operate_binary(k, k), neg_one);

    // ij = k
    assert_eq!(mul.operate_binary(i, j), k);
    assert_eq!(mul.operate_binary(j, i), Quaternion { re: 0.0, i: 0.0, j: 0.0, k: -1.0 });

    // jk = i
    assert_eq!(mul.operate_binary(j, k), i);
    assert_eq!(mul.operate_binary(k, j), Quaternion { re: 0.0, i: -1.0, j: 0.0, k: 0.0 });

    // ki = j
    assert_eq!(mul.operate_binary(k, i), j);
    assert_eq!(mul.operate_binary(i, k), Quaternion { re: 0.0, i: 0.0, j: -1.0, k: 0.0 });
}

#[test]
fn test_quaternion_multiplication_associativity() {
    let mul = QuaternionMultiplication;

    let a = Quaternion { re: 1.0, i: 2.0, j: 3.0, k: 4.0 };
    let b = Quaternion { re: -1.0, i: 1.0, j: -1.0, k: 1.0 };
    let c = Quaternion { re: 0.5, i: -2.0, j: 0.0, k: 3.0 };

    let ab = mul.operate_binary(a, b);
    let bc = mul.operate_binary(b, c);

    let ab_c = mul.operate_binary(ab, c);
    let a_bc = mul.operate_binary(a, bc);

    // Allow some floating-point tolerance
    fn approx_eq(q1: Quaternion, q2: Quaternion) -> bool {
        (q1.re - q2.re).abs() < 1e-6 &&
        (q1.i - q2.i).abs() < 1e-6 &&
        (q1.j - q2.j).abs() < 1e-6 &&
        (q1.k - q2.k).abs() < 1e-6
    }

    assert!(
        approx_eq(ab_c, a_bc),
        "Associativity failed: (ab)c = {:?}, a(bc) = {:?}",
        ab_c, a_bc
    );
}
