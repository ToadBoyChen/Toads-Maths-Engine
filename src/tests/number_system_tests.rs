// mod libs {
//     include!("../lib/algebras.rs");
//     include!("../lib/specialised_algebras.rs");
// }

// use crate::libs::*;

macro_rules! assert_complex_eq {
    ($a:expr, $b:expr) => {
        assert!(($a.re - $b.re).abs() < 1e-6 && ($a.im - $b.im).abs() < 1e-6,
      
            "Complex mismatch: {:?} != {:?}", $a, $b);
    };
} 

macro_rules! assert_quaternion_eq {
    ($a:expr, $b:expr) => {
        assert!(($a.re - $b.re).abs() < 1e-6 && ($a.i - $b.i).abs() < 1e-6 &&
                ($a.j - $b.j).abs() < 1e-6 && ($a.k - $b.k).abs() < 1e-6,
            "Quaternion mismatch: {:?} != {:?}", $a, $b);
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
fn test_quaternion_addition_group_laws() {
    let add = QuaternionAddition;
    let mut rng = rand::thread_rng();

    for _ in 0..100 {
        let a = Quaternion { re: rng.gen_range(-100.0..100.0), i: rng.gen_range(-100.0..100.0), j: rng.gen_range(-100.0..100.0), k: rng.gen_range(-100.0..100.0) };
        let b = Quaternion { re: rng.gen_range(-100.0..100.0), i: rng.gen_range(-100.0..100.0), j: rng.gen_range(-100.0..100.0), k: rng.gen_range(-100.0..100.0) };
        let c = Quaternion { re: rng.gen_range(-100.0..100.0), i: rng.gen_range(-100.0..100.0), j: rng.gen_range(-100.0..100.0), k: rng.gen_range(-100.0..100.0) };

        // Associativity
        assert_quaternion_eq!(
            add.operate_binary(add.operate_binary(a, b), c),
            add.operate_binary(a, add.operate_binary(b, c))
        );

        // Commutativity
        assert_quaternion_eq!(add.operate_binary(a, b), add.operate_binary(b, a));

        // Identity
        assert_quaternion_eq!(add.operate_binary(a, Quaternion { re: 0.0, i: 0.0, j: 0.0, k: 0.0 }), a);
        assert_quaternion_eq!(add.operate_binary(Quaternion { re: 0.0, i: 0.0, j: 0.0, k: 0.0 }, a), a);

        // Inverse
        let inv = Quaternion { re: -a.re, i: -a.i, j: -a.j, k: -a.k };
        assert_quaternion_eq!(add.operate_binary(a, inv), Quaternion { re: 0.0, i: 0.0, j: 0.0, k: 0.0 });
    }

    let a = Quaternion { re: 5.0, i: 0.0, j: 7.0, k: 0.0 };
    assert_quaternion_eq!(add.operate(&[a]), a);
}


#[test]
fn test_quaternion_multiplication_group_laws() {
    let mul = QuaternionMultiplication;
    let mut rng = rand::thread_rng();
    let one = Quaternion { re: 1.0, i: 0.0, j: 0.0, k: 0.0 };

    for _ in 0..100 {
        let a = Quaternion { re: rng.gen_range(-100.0..100.0), i: rng.gen_range(-100.0..100.0), j: rng.gen_range(-100.0..100.0), k: rng.gen_range(-100.0..100.0) };
        let b = Quaternion { re: rng.gen_range(-100.0..100.0), i: rng.gen_range(-100.0..100.0), j: rng.gen_range(-100.0..100.0), k: rng.gen_range(-100.0..100.0) };
        let c = Quaternion { re: rng.gen_range(-100.0..100.0), i: rng.gen_range(-100.0..100.0), j: rng.gen_range(-100.0..100.0), k: rng.gen_range(-100.0..100.0) };

        // Associativity
        assert_quaternion_eq!(
            mul.operate_binary(mul.operate_binary(a, b), c),
            mul.operate_binary(a, mul.operate_binary(b, c))
        );

        // Commutativity - ALL quaternion multiplication is not commutative
        assert_ne!(mul.operate_binary(a, b), mul.operate_binary(b, a), "Quaternions are unexpectedly commutative: {:?} * {:?} == {:?} * {:?}", a, b, b, a);

        // Identity
        assert_quaternion_eq!(mul.operate_binary(a, one), a);
        assert_quaternion_eq!(mul.operate_binary(one, a), a);

        // Inverse
        let norm_sq = a.re * a.re + a.i * a.i + a.j * a.j + a.k * a.k;
        let inv = Quaternion {
            re: a.re / norm_sq,
            i: -a.i / norm_sq,
            j: -a.j / norm_sq,
            k: -a.k / norm_sq,
        };
        assert_quaternion_eq!(mul.operate_binary(a, inv), one);
    }

    assert_quaternion_eq!(mul.operate(&[]), one);
    let a = Quaternion { re: 5.0, i: 0.0, j: 7.0, k: 0.0 };
    assert_quaternion_eq!(mul.operate(&[a]), a);
}


// Complex addition: use approx_eq
#[test]
fn test_complex_addition_group_laws() {
    let add = ComplexAddition;
    let mut rng = rand::thread_rng();

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
        assert_complex_eq!(add.operate_binary(a, Complex { re: 0.0, im: 0.0 }), a);
        assert_complex_eq!(add.operate_binary(Complex { re: 0.0, im: 0.0 }, a), a);

        // Inverse
        let inv = Complex { re: -a.re, im: -a.im };
        assert_complex_eq!(add.operate_binary(a, inv), Complex { re: 0.0, im: 0.0 });
    }

    let a = Complex { re: 5.0, im: 7.0 };
    assert_complex_eq!(add.operate(&[a]), a);
}


#[test]
fn test_complex_multiplication_group_laws() {
    let mul = ComplexMultiplication;
    let mut rng = rand::thread_rng();
    let one = Complex { re: 1.0, im: 0.0 };

    for _ in 0..100 {
        let a = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        let b = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        let c = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };

        // Associativity
        assert_complex_eq!(
            mul.operate_binary(mul.operate_binary(a, b), c),
            mul.operate_binary(a, mul.operate_binary(b, c))
        );

        // Commutativity
        assert_complex_eq!(mul.operate_binary(a, b), mul.operate_binary(b, a));

        // Identity
        assert_complex_eq!(mul.operate_binary(a, one), a);
        assert_complex_eq!(mul.operate_binary(one, a), a);

        // Inverse
        let norm_sq = a.re * a.re + a.im * a.im;
        let inv = Complex {
            re: a.re / norm_sq,
            im: -a.im / norm_sq,
        };
        assert_complex_eq!(mul.operate_binary(a, inv), one);
    }

    assert_complex_eq!(mul.operate(&[]), one);
    let a = Complex { re: 5.0, im: 7.0 };
    assert_complex_eq!(mul.operate(&[a]), a);
}

#[test]
fn test_complex_field_laws() {
    let add = ComplexAddition;
    let mul = ComplexMultiplication;
    let mut rng = rand::thread_rng();
    let one = Complex { re: 1.0, im: 0.0 };
    let zero = Complex { re: 0.0, im: 0.0 };

    for _ in 0..100 {
        let a = Complex { re: rng.gen_range(-50.0..50.0), im: rng.gen_range(-50.0..50.0) };
        let b = Complex { re: rng.gen_range(-50.0..50.0), im: rng.gen_range(-50.0..50.0) };
        let c = Complex { re: rng.gen_range(-50.0..50.0), im: rng.gen_range(-50.0..50.0) };

        // Distributivity
        let lhs = mul.operate_binary(a, add.operate_binary(b, c));
        let rhs = add.operate_binary(mul.operate_binary(a, b), mul.operate_binary(a, c));
        assert_complex_eq!(lhs, rhs);

        // Multiplicative inverse exists (if non-zero)
        if a.re != 0.0 || a.im != 0.0 {
            let norm_sq = a.re * a.re + a.im * a.im;
            let inv = Complex {
                re: a.re / norm_sq,
                im: -a.im / norm_sq,
            };
            let prod = mul.operate_binary(a, inv);
            assert_complex_eq!(prod, one);
        }
    }
}

#[test]
fn test_octonion_structures() {
    let mut rng = rand::thread_rng();

    // Test Octonion Addition
    let add = OctonionAddition;
    let a = Octonion {
        re: rng.gen_range(-100.0..100.0),
        i: rng.gen_range(-100.0..100.0),
        j: rng.gen_range(-100.0..100.0),
        k: rng.gen_range(-100.0..100.0),
        l: rng.gen_range(-100.0..100.0),
        m: rng.gen_range(-100.0..100.0),
        n: rng.gen_range(-100.0..100.0),
        o: rng.gen_range(-100.0..100.0),
    };
    let b = Octonion {
        re: rng.gen_range(-100.0..100.0),
        i: rng.gen_range(-100.0..100.0),
        j: rng.gen_range(-100.0..100.0),
        k: rng.gen_range(-100.0..100.0),
        l: rng.gen_range(-100.0..100.0),
        m: rng.gen_range(-100.0..100.0),
        n: rng.gen_range(-100.0..100.0),
        o: rng.gen_range(-100.0..100.0),
    };

    let sum = add.operate_binary(a, b);
    assert_eq!(
        sum,
        Octonion {
            re: a.re + b.re,
            i: a.i + b.i,
            j: a.j + b.j,
            k: a.k + b.k,
            l: a.l + b.l,
            m: a.m + b.m,
            n: a.n + b.n,
            o: a.o + b.o,
        }
    );

    // Test Octonion Multiplication
    let mul = OctonionMultiplication;
    let product = mul.operate_binary(a, b);
    assert_ne!(product, mul.operate_binary(b, a), "Octonion multiplication is unexpectedly commutative");

    // Test Quaternion Multiplication
    let quaternion_mul = QuaternionMultiplication;
    let q1 = Quaternion {
        re: rng.gen_range(-100.0..100.0),
        i: rng.gen_range(-100.0..100.0),
        j: rng.gen_range(-100.0..100.0),
        k: rng.gen_range(-100.0..100.0),
    };
    let q2 = Quaternion {
        re: rng.gen_range(-100.0..100.0),
        i: rng.gen_range(-100.0..100.0),
        j: rng.gen_range(-100.0..100.0),
        k: rng.gen_range(-100.0..100.0),
    };

    let quaternion_product = quaternion_mul.operate_binary(q1, q2);
    assert_ne!(quaternion_product, quaternion_mul.operate_binary(q2, q1), "Quaternion multiplication is unexpectedly commutative");
}

#[test]
fn test_non_commutativity_quaternions() {
    let mul = QuaternionMultiplication;
    let mut rng = rand::thread_rng();

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

    let ab = mul.operate_binary(a, b);
    let ba = mul.operate_binary(b, a);

    assert_ne!(ab, ba, "Quaternion multiplication is unexpectedly commutative");
}

#[test]
#[should_panic(expected = "Cannot invert a zero quaternion")]
fn test_zero_quaternion_inverse() {
    let mul = QuaternionMultiplication;

    // Define the zero quaternion
    let zero_quaternion = Quaternion {
        re: 0.0,
        i: 0.0,
        j: 0.0,
        k: 0.0,
    };

    // Attempt to compute the inverse, which should panic
    mul.inverse(zero_quaternion);
}