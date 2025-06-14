macro_rules! assert_complex_eq {
    ($a:expr, $b:expr, $tolerance:expr) => {
        assert!(($a.re - $b.re).abs() < $tolerance && ($a.im - $b.im).abs() < $tolerance,
            "Complex mismatch: {:?} != {:?} (tolerance: {})", $a, $b, $tolerance);
    };
    ($a:expr, $b:expr) => {
        assert_complex_eq!($a, $b, 1e-10);
    };
} 

macro_rules! assert_quaternion_eq {
    ($a:expr, $b:expr, $tolerance:expr) => {
        assert!(($a.re - $b.re).abs() < $tolerance && ($a.i - $b.i).abs() < $tolerance &&
                ($a.j - $b.j).abs() < $tolerance && ($a.k - $b.k).abs() < $tolerance,
            "Quaternion mismatch: {:?} != {:?} (tolerance: {})", $a, $b, $tolerance);
    };
    ($a:expr, $b:expr) => {
        assert_quaternion_eq!($a, $b, 1e-10);
    };
}

macro_rules! assert_octonion_eq {
    ($a:expr, $b:expr, $tolerance:expr) => {
        assert!(($a.re - $b.re).abs() < $tolerance && ($a.i - $b.i).abs() < $tolerance &&
                ($a.j - $b.j).abs() < $tolerance && ($a.k - $b.k).abs() < $tolerance &&
                ($a.l - $b.l).abs() < $tolerance && ($a.m - $b.m).abs() < $tolerance &&
                ($a.n - $b.n).abs() < $tolerance && ($a.o - $b.o).abs() < $tolerance,
            "Octonion mismatch: {:?} != {:?} (tolerance: {})", $a, $b, $tolerance);
    };
    ($a:expr, $b:expr) => {
        assert_octonion_eq!($a, $b, 1e-10);
    };
}

// ========== INTEGER ADDITION GROUP TESTS ==========

#[test]
fn test_integer_addition_group_axioms() {
    let add = IntegerAddition;
    let test_values = [-1000, -100, -1, 0, 1, 100, 1000];

    // Test closure (with wrapping arithmetic)
    for &a in &test_values {
        for &b in &test_values {
            let result = add.operate_binary(a, b);
            assert_eq!(result, a.wrapping_add(b), "Closure failed for {} + {}", a, b);
        }
    }

    // Test associativity: (a + b) + c = a + (b + c)
    for &a in &test_values {
        for &b in &test_values {
            for &c in &test_values {
                let left = add.operate_binary(add.operate_binary(a, b), c);
                let right = add.operate_binary(a, add.operate_binary(b, c));
                assert_eq!(left, right, "Associativity failed for ({} + {}) + {} vs {} + ({} + {})", a, b, c, a, b, c);
            }
        }
    }

    // Test commutativity: a + b = b + a (abelian group)
    for &a in &test_values {
        for &b in &test_values {
            assert_eq!(add.operate_binary(a, b), add.operate_binary(b, a), 
                "Commutativity failed for {} + {} vs {} + {}", a, b, b, a);
        }
    }

    // Test identity element: a + 0 = a
    assert_eq!(add.identity(), 0);
    for &a in &test_values {
        assert_eq!(add.operate_binary(a, 0), a, "Right identity failed for {}", a);
        assert_eq!(add.operate_binary(0, a), a, "Left identity failed for {}", a);
    }

    // Test inverse element: a + (-a) = 0
    for &a in &test_values {
        let inverse = add.inverse(a);
        assert_eq!(inverse, -a, "Inverse calculation failed for {}", a);
        assert_eq!(add.operate_binary(a, inverse), 0, "Inverse property failed for {}", a);
    }

    // Test operation on collections
    assert_eq!(add.operate(&[]), 0, "Empty collection should return identity");
    assert_eq!(add.operate(&[42]), 42, "Single element should return itself");
    assert_eq!(add.operate(&[1, 2, 3, 4]), 10, "Sum of [1,2,3,4] should be 10");
}

#[test]
fn test_integer_addition_edge_cases() {
    let add = IntegerAddition;

    // Test overflow behavior
    assert_eq!(add.operate_binary(i32::MAX, 1), i32::MIN, "Overflow should wrap");
    assert_eq!(add.operate_binary(i32::MIN, -1), i32::MAX, "Underflow should wrap");
    
    // Test large sequences
    let large_seq: Vec<i32> = (1..=1000).collect();
    let expected_sum = (1000 * 1001) / 2; // Gauss formula
    assert_eq!(add.operate(&large_seq), expected_sum, "Large sequence sum failed");
}

// ========== INTEGER MULTIPLICATION MONOID TESTS ==========

#[test]
fn test_integer_multiplication_monoid_axioms() {
    let mul = IntegerMultiplication;
    let test_values = [-10, -2, -1, 0, 1, 2, 10];

    // Test closure (with wrapping arithmetic)
    for &a in &test_values {
        for &b in &test_values {
            let result = mul.operate_binary(a, b);
            assert_eq!(result, a.wrapping_mul(b), "Closure failed for {} * {}", a, b);
        }
    }

    // Test associativity: (a * b) * c = a * (b * c)
    for &a in &test_values {
        for &b in &test_values {
            for &c in &test_values {
                let left = mul.operate_binary(mul.operate_binary(a, b), c);
                let right = mul.operate_binary(a, mul.operate_binary(b, c));
                assert_eq!(left, right, "Associativity failed for ({} * {}) * {} vs {} * ({} * {})", a, b, c, a, b, c);
            }
        }
    }

    // Test commutativity: a * b = b * a (abelian monoid)
    for &a in &test_values {
        for &b in &test_values {
            assert_eq!(mul.operate_binary(a, b), mul.operate_binary(b, a),
                "Commutativity failed for {} * {} vs {} * {}", a, b, b, a);
        }
    }

    // Test identity element: a * 1 = a
    for &a in &test_values {
        assert_eq!(mul.operate_binary(a, 1), a, "Right identity failed for {}", a);
        assert_eq!(mul.operate_binary(1, a), a, "Left identity failed for {}", a);
    }

    // Test operation on collections
    assert_eq!(mul.operate(&[]), 1, "Empty collection should return identity");
    assert_eq!(mul.operate(&[7]), 7, "Single element should return itself");
    assert_eq!(mul.operate(&[2, 3, 4]), 24, "Product of [2,3,4] should be 24");
}

#[test]
fn test_integer_multiplication_special_cases() {
    let mul = IntegerMultiplication;

    // Test zero element (absorbing element)
    for &a in &[-100, -1, 0, 1, 100] {
        assert_eq!(mul.operate_binary(a, 0), 0, "Zero should absorb: {} * 0", a);
        assert_eq!(mul.operate_binary(0, a), 0, "Zero should absorb: 0 * {}", a);
    }

    // Test negative numbers
    assert_eq!(mul.operate_binary(-2, -3), 6, "Negative * negative should be positive");
    assert_eq!(mul.operate_binary(-2, 3), -6, "Negative * positive should be negative");
    assert_eq!(mul.operate_binary(2, -3), -6, "Positive * negative should be negative");
}

// ========== COMPLEX NUMBER TESTS ==========

#[test]
fn test_complex_addition_group_axioms() {
    let add = ComplexAddition;
    let mut rng = rand::thread_rng();

    // Test with specific values
    let test_values = [
        Complex { re: 0.0, im: 0.0 },
        Complex { re: 1.0, im: 0.0 },
        Complex { re: 0.0, im: 1.0 },
        Complex { re: 1.0, im: 1.0 },
        Complex { re: -1.0, im: -1.0 },
        Complex { re: 3.0, im: -4.0 },
    ];

    for &a in &test_values {
        for &b in &test_values {
            for &c in &test_values {
                // Test associativity: (a + b) + c = a + (b + c)
                let left = add.operate_binary(add.operate_binary(a, b), c);
                let right = add.operate_binary(a, add.operate_binary(b, c));
                assert_complex_eq!(left, right);

                // Test commutativity: a + b = b + a
                assert_complex_eq!(add.operate_binary(a, b), add.operate_binary(b, a));
            }

            // Test identity: a + 0 = a
            let zero = Complex { re: 0.0, im: 0.0 };
            assert_complex_eq!(add.operate_binary(a, zero), a);
            assert_complex_eq!(add.operate_binary(zero, a), a);

            // Test inverse: a + (-a) = 0
            let inverse = add.inverse(a);
            assert_complex_eq!(add.operate_binary(a, inverse), zero);
        }
    }

    // Test with random values for robustness
    for _ in 0..50 {
        let a = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        let b = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        let c = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };

        // Test associativity
        let left = add.operate_binary(add.operate_binary(a, b), c);
        let right = add.operate_binary(a, add.operate_binary(b, c));
        assert_complex_eq!(left, right);

        // Test commutativity
        assert_complex_eq!(add.operate_binary(a, b), add.operate_binary(b, a));
    }
}

#[test]
fn test_complex_multiplication_group_axioms() {
    let mul = ComplexMultiplication;
    let mut rng = rand::thread_rng();

    // Test with specific values
    let test_values = [
        Complex { re: 1.0, im: 0.0 },   // Real unit
        Complex { re: 0.0, im: 1.0 },   // Imaginary unit
        Complex { re: -1.0, im: 0.0 },  // Negative real
        Complex { re: 0.0, im: -1.0 },  // Negative imaginary
        Complex { re: 1.0, im: 1.0 },   // First quadrant
        Complex { re: 3.0, im: 4.0 },   // Pythagorean triple
    ];

    for &a in &test_values {
        for &b in &test_values {
            for &c in &test_values {
                // Test associativity: (a * b) * c = a * (b * c)
                let left = mul.operate_binary(mul.operate_binary(a, b), c);
                let right = mul.operate_binary(a, mul.operate_binary(b, c));
                assert_complex_eq!(left, right);

                // Test commutativity: a * b = b * a
                assert_complex_eq!(mul.operate_binary(a, b), mul.operate_binary(b, a));
            }

            // Test identity: a * 1 = a
            let one = Complex { re: 1.0, im: 0.0 };
            assert_complex_eq!(mul.operate_binary(a, one), a);
            assert_complex_eq!(mul.operate_binary(one, a), a);

            // Test inverse: a * a⁻¹ = 1 (for non-zero a)
            if a.re != 0.0 || a.im != 0.0 {
                let inverse = mul.inverse(a);
                let product = mul.operate_binary(a, inverse);
                assert_complex_eq!(product, one);
            }
        }
    }

    // Test random values
    for _ in 0..50 {
        let a = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        let b = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        
        // Skip if too close to zero
        if (a.re * a.re + a.im * a.im) < 1e-10 || (b.re * b.re + b.im * b.im) < 1e-10 {
            continue;
        }

        // Test commutativity
        assert_complex_eq!(mul.operate_binary(a, b), mul.operate_binary(b, a));

        // Test inverse
        let inv_a = mul.inverse(a);
        let product = mul.operate_binary(a, inv_a);
        assert_complex_eq!(product, Complex { re: 1.0, im: 0.0 });
    }
}

#[test]
fn test_complex_field_properties() {
    let add = ComplexAddition;
    let mul = ComplexMultiplication;
    let mut rng = rand::thread_rng();

    for _ in 0..50 {
        let a = Complex { re: rng.gen_range(-50.0..50.0), im: rng.gen_range(-50.0..50.0) };
        let b = Complex { re: rng.gen_range(-50.0..50.0), im: rng.gen_range(-50.0..50.0) };
        let c = Complex { re: rng.gen_range(-50.0..50.0), im: rng.gen_range(-50.0..50.0) };

        // Test distributivity: a * (b + c) = a * b + a * c
        let left = mul.operate_binary(a, add.operate_binary(b, c));
        let right = add.operate_binary(mul.operate_binary(a, b), mul.operate_binary(a, c));
        assert_complex_eq!(left, right);

        // Test right distributivity: (a + b) * c = a * c + b * c
        let left_right = mul.operate_binary(add.operate_binary(a, b), c);
        let right_right = add.operate_binary(mul.operate_binary(a, c), mul.operate_binary(b, c));
        assert_complex_eq!(left_right, right_right);
    }
}

#[test]
fn test_complex_edge_cases() {
    let add = ComplexAddition;
    let mul = ComplexMultiplication;

    // Test zero element
    let zero = Complex { re: 0.0, im: 0.0 };
    assert_complex_eq!(add.operate_binary(zero, zero), zero);
    assert_complex_eq!(mul.operate_binary(zero, zero), zero);

    // Test identity element
    let one = Complex { re: 1.0, im: 0.0 };
    assert_complex_eq!(add.operate_binary(one, zero), one);
    assert_complex_eq!(mul.operate_binary(one, one), one);

    // Test inverse element
    let neg_one = Complex { re: -1.0, im: 0.0 };
    assert_complex_eq!(add.operate_binary(one, neg_one), zero);
}
#[test]
fn test_complex_randomized_operations() {
    let add = ComplexAddition;
    let mul = ComplexMultiplication;
    let mut rng = rand::thread_rng();

    for _ in 0..100 {
        let a = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };
        let b = Complex { re: rng.gen_range(-100.0..100.0), im: rng.gen_range(-100.0..100.0) };

        // Test addition
        let sum = add.operate_binary(a, b);
        assert_complex_eq!(sum, Complex { re: a.re + b.re, im: a.im + b.im });

        // Test multiplication
        let product = mul.operate_binary(a, b);
        assert_complex_eq!(product, Complex {
            re: a.re * b.re - a.im * b.im,
            im: a.re * b.im + a.im * b.re,
        });
    }
}

#[test]
#[should_panic(expected = "division by zero")] // Or whatever the expected panic message is
fn test_complex_inverse_of_zero_panics() {
    let mul = ComplexMultiplication;
    let zero = Complex { re: 0.0, im: 0.0 };
    mul.inverse(zero); // This line should cause a panic
}