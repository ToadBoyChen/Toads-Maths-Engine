mod libs {
    include!("../lib/algebras.rs");
    include!("../lib/specialised_algebras.rs");
}

use crate::libs::*;
use rand::Rng;

// ========== EUCLIDEAN SPACE TESTS ==========

#[test]
fn test_euclidean_space_vector_properties() {
    let space = EuclideanSpace::new(3);

    // Standard orthonormal basis vectors
    let e1 = space.vector(&[1.0, 0.0, 0.0]);
    let e2 = space.vector(&[0.0, 1.0, 0.0]);
    let e3 = space.vector(&[0.0, 0.0, 1.0]);

    // Test orthogonality of basis vectors
    assert_eq!(space.dot(&e1, &e2), 0.0, "e1 and e2 should be orthogonal");
    assert_eq!(space.dot(&e1, &e3), 0.0, "e1 and e3 should be orthogonal");
    assert_eq!(space.dot(&e2, &e3), 0.0, "e2 and e3 should be orthogonal");

    // Test unit lengths
    assert!((space.norm(&e1) - 1.0).abs() < 1e-10, "e1 should have unit length");
    assert!((space.norm(&e2) - 1.0).abs() < 1e-10, "e2 should have unit length");
    assert!((space.norm(&e3) - 1.0).abs() < 1e-10, "e3 should have unit length");

    // Test Pythagorean theorem
    let v = space.vector(&[3.0, 4.0, 0.0]);
    assert!((space.norm(&v) - 5.0).abs() < 1e-10, "3-4-5 triangle should have norm 5");

    // Test vector addition commutativity
    let u = space.vector(&[1.0, 2.0, 3.0]);
    let w = space.vector(&[4.0, -1.0, 2.0]);
    assert_eq!(space.add(&u, &w), space.add(&w, &u), "Vector addition should be commutative");

    // Test scalar multiplication distributivity: a(u + v) = au + av
    let scalar = 2.5;
    let sum = space.add(&u, &w);
    let scaled_sum = space.scalar_mul(&sum, scalar);
    let sum_of_scaled = space.add(&space.scalar_mul(&u, scalar), &space.scalar_mul(&w, scalar));
    
    for (a, b) in scaled_sum.iter().zip(sum_of_scaled.iter()) {
        assert!((a - b).abs() < 1e-10, "Scalar multiplication should distribute over addition");
    }

    // Test dot product bilinearity: (au + bv) · w = a(u·w) + b(v·w)
    let a = 2.0;
    let b = 3.0;
    let au_plus_bv = space.add(&space.scalar_mul(&u, a), &space.scalar_mul(&w, b));
    let left_side = space.dot(&au_plus_bv, &e1);
    let right_side = a * space.dot(&u, &e1) + b * space.dot(&w, &e1);
    assert!((left_side - right_side).abs() < 1e-10, "Dot product should be bilinear");
}

#[test]
fn test_euclidean_space_geometric_properties() {
    let space = EuclideanSpace::new(2);

    // Test Cauchy-Schwarz inequality: |u·v| ≤ ||u|| ||v||
    let u = space.vector(&[3.0, 4.0]);
    let v = space.vector(&[1.0, 2.0]);
    
    let dot_product = space.dot(&u, &v);
    let norm_u = space.norm(&u);
    let norm_v = space.norm(&v);
    
    assert!(dot_product.abs() <= norm_u * norm_v + 1e-10, 
        "Cauchy-Schwarz inequality should hold: |{}| ≤ {} * {}", dot_product, norm_u, norm_v);

    // Test parallelogram law: ||u + v||² + ||u - v||² = 2(||u||² + ||v||²)
    let sum = space.add(&u, &v);
    let diff = space.add(&u, &space.scalar_mul(&v, -1.0));
    
    let left_side = space.norm(&sum).powi(2) + space.norm(&diff).powi(2);
    let right_side = 2.0 * (space.norm(&u).powi(2) + space.norm(&v).powi(2));
    
    assert!((left_side - right_side).abs() < 1e-10, 
        "Parallelogram law should hold: {} ≈ {}", left_side, right_side);
}

// ========== ELEMENT-WISE VECTOR OPERATIONS TESTS ==========

#[test]
fn test_element_wise_addition_group_properties() {
    let add = NMatrixAddition::zero();
    let mut rng = rand::thread_rng();

    for dimension in 1..=10 {
        // Create test vectors with random values
        let a = NMatrixAddition { 
            elements: (0..dimension).map(|_| rng.gen_range(-100.0..100.0)).collect() 
        };
        let b = NMatrixAddition { 
            elements: (0..dimension).map(|_| rng.gen_range(-100.0..100.0)).collect() 
        };
        let c = NMatrixAddition { 
            elements: (0..dimension).map(|_| rng.gen_range(-100.0..100.0)).collect() 
        };

        // Test closure: result has same dimension
        let sum = add.operate_binary(a.clone(), b.clone());
        assert_eq!(sum.elements.len(), dimension, "Addition should preserve dimension");

        // Test associativity: (a + b) + c = a + (b + c)
        let left_assoc = add.operate_binary(add.operate_binary(a.clone(), b.clone()), c.clone());
        let right_assoc = add.operate_binary(a.clone(), add.operate_binary(b.clone(), c.clone()));
        
        for (l, r) in left_assoc.elements.iter().zip(right_assoc.elements.iter()) {
            assert!((l - r).abs() < 1e-10, "Addition should be associative");
        }

        // Test commutativity: a + b = b + a
        let ab = add.operate_binary(a.clone(), b.clone());
        let ba = add.operate_binary(b.clone(), a.clone());
        assert_eq!(ab, ba, "Addition should be commutative");

        // Test identity: a + 0 = a
        let zero = NMatrixAddition { elements: vec![0.0; dimension] };
        let identity_result = add.operate_binary(a.clone(), zero.clone());
        assert_eq!(identity_result, a, "Zero should be additive identity");

        // Test inverse: a + (-a) = 0
        let inverse = add.inverse(a.clone());
        let inverse_sum = add.operate_binary(a.clone(), inverse);
        
        for &val in &inverse_sum.elements {
            assert!(val.abs() < 1e-10, "Inverse should sum to zero");
        }
    }
}

#[test]
fn test_element_wise_multiplication_group_properties() {
    let mul = NMatrixMultiplication::zero();
    let mut rng = rand::thread_rng();

    for dimension in 1..=10 {
        // Create test vectors with non-zero random values
        let a = NMatrixMultiplication { 
            elements: (0..dimension).map(|_| {
                let val: f64 = rng.gen_range(-100.0..100.0);
                if val.abs() < 1e-6 { 1.0 } else { val } // Avoid near-zero values
            }).collect() 
        };
        let b = NMatrixMultiplication { 
            elements: (0..dimension).map(|_| {
                let val: f64 = rng.gen_range(-100.0..100.0);
                if val.abs() < 1e-6 { 1.0 } else { val }
            }).collect() 
        };
        let c = NMatrixMultiplication { 
            elements: (0..dimension).map(|_| {
                let val: f64 = rng.gen_range(-100.0..100.0);
                if val.abs() < 1e-6 { 1.0 } else { val }
            }).collect() 
        };

        // Test closure: result has same dimension
        let product = mul.operate_binary(a.clone(), b.clone());
        assert_eq!(product.elements.len(), dimension, "Multiplication should preserve dimension");

        // Test associativity: (a * b) * c = a * (b * c)
        let left_assoc = mul.operate_binary(mul.operate_binary(a.clone(), b.clone()), c.clone());
        let right_assoc = mul.operate_binary(a.clone(), mul.operate_binary(b.clone(), c.clone()));
        
        for (l, r) in left_assoc.elements.iter().zip(right_assoc.elements.iter()) {
            assert!((l - r).abs() < 1e-10, "Multiplication should be associative");
        }

        // Test commutativity: a * b = b * a
        let ab = mul.operate_binary(a.clone(), b.clone());
        let ba = mul.operate_binary(b.clone(), a.clone());
        assert_eq!(ab, ba, "Element-wise multiplication should be commutative");

        // Test identity: a * 1 = a
        let ones = NMatrixMultiplication { elements: vec![1.0; dimension] };
        let identity_result = mul.operate_binary(a.clone(), ones);
        assert_eq!(identity_result, a, "Ones vector should be multiplicative identity");

        // Test inverse: a * a⁻¹ = 1
        let inverse = mul.inverse(a.clone());
        let inverse_product = mul.operate_binary(a.clone(), inverse);
        
        for &val in &inverse_product.elements {
            assert!((val - 1.0).abs() < 1e-10, "Inverse should multiply to one: got {}", val);
        }
    }
}

#[test]
fn test_vector_operations_distributivity() {
    let add = NMatrixAddition::zero();
    let mut rng = rand::thread_rng();

    // Test distributivity of scalar multiplication over addition
    for _ in 0..20 {
        let dimension = rng.gen_range(2..20);
        let scalar = rng.gen_range(-10.0..10.0);
        
        let a = NMatrixAddition { 
            elements: (0..dimension).map(|_| rng.gen_range(-100.0..100.0)).collect() 
        };
        let b = NMatrixAddition { 
            elements: (0..dimension).map(|_| rng.gen_range(-100.0..100.0)).collect() 
        };

        // Test: scalar * (a + b) = scalar * a + scalar * b
        let sum_ab = add.operate_binary(a.clone(), b.clone());
        let scalar_sum: Vec<f64> = sum_ab.elements.iter().map(|&x| scalar * x).collect();
        
        let scalar_a: Vec<f64> = a.elements.iter().map(|&x| scalar * x).collect();
        let scalar_b: Vec<f64> = b.elements.iter().map(|&x| scalar * x).collect();
        let sum_scalars: Vec<f64> = scalar_a.iter().zip(scalar_b.iter()).map(|(&x, &y)| x + y).collect();

        for (left, right) in scalar_sum.iter().zip(sum_scalars.iter()) {
            assert!((left - right).abs() < 1e-10, 
                "Scalar multiplication should distribute over addition: {} ≠ {}", left, right);
        }
    }
}

#[test]
fn test_edge_cases_and_special_values() {
    let add = NMatrixAddition::zero();
    let mul = NMatrixMultiplication::zero();

    // Test empty vectors
    let empty1 = NMatrixAddition { elements: vec![] };
    let empty2 = NMatrixAddition { elements: vec![] };
    let empty_sum = add.operate_binary(empty1, empty2);
    assert_eq!(empty_sum.elements.len(), 0, "Empty vector operations should work");

    // Test single element vectors
    let single1 = NMatrixAddition { elements: vec![5.0] };
    let single2 = NMatrixAddition { elements: vec![3.0] };
    let single_sum = add.operate_binary(single1, single2);
    assert_eq!(single_sum.elements, vec![8.0], "Single element addition should work");

    // Test very large vectors
    let large_size = 10000;
    let large1 = NMatrixAddition { elements: vec![1.0; large_size] };
    let large2 = NMatrixAddition { elements: vec![2.0; large_size] };
    let large_sum = add.operate_binary(large1, large2);
    assert_eq!(large_sum.elements.len(), large_size, "Large vector operations should work");
    assert!(large_sum.elements.iter().all(|&x| x == 3.0), "Large vector sum should be correct");

    // Test precision with very small numbers
    let small1 = NMatrixAddition { elements: vec![1e-15, 2e-15] };
    let small2 = NMatrixAddition { elements: vec![3e-15, 4e-15] };
    let small_sum = add.operate_binary(small1, small2);
    assert!((small_sum.elements[0] - 4e-15).abs() < 1e-20, "Small number precision should be maintained");
    assert!((small_sum.elements[1] - 6e-15).abs() < 1e-20, "Small number precision should be maintained");
}

#[test]
fn test_operation_chaining() {
    let add = NMatrixAddition::zero();
    let mul = NMatrixMultiplication::zero();

    // Test chaining multiple additions
    let vectors = vec![
        NMatrixAddition { elements: vec![1.0, 2.0, 3.0] },
        NMatrixAddition { elements: vec![4.0, 5.0, 6.0] },
        NMatrixAddition { elements: vec![7.0, 8.0, 9.0] },
        NMatrixAddition { elements: vec![10.0, 11.0, 12.0] },
    ];

    let sum_result = add.operate(&vectors);
    assert_eq!(sum_result.elements, vec![22.0, 26.0, 30.0], "Chained addition should work");

    // Test chaining multiple multiplications
    let mul_vectors = vec![
        NMatrixMultiplication { elements: vec![2.0, 3.0] },
        NMatrixMultiplication { elements: vec![4.0, 2.0] },
        NMatrixMultiplication { elements: vec![0.5, 2.0] },
    ];

    let mul_result = mul.operate(&mul_vectors);
    assert_eq!(mul_result.elements, vec![4.0, 12.0], "Chained multiplication should work");
}