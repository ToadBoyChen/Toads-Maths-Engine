mod libs {
    include!("../lib/algebras.rs");
    include!("../lib/specialised_algebras.rs");
}

use crate::libs::*;
use rand::Rng;


#[test]
fn test_euclidean_space_properties() {
    let space = EuclideanSpace::new(3); // 3D space

    let v1 = space.vector(&[1.0, 0.0, 0.0]);
    let v2 = space.vector(&[0.0, 1.0, 0.0]);

    // Dot product
    assert_eq!(space.dot(&v1, &v2), 0.0);

    // Norm
    assert!((space.norm(&v1) - 1.0).abs() < 1e-6);

    // Vector addition
    let sum = space.add(&v1, &v2);
    assert_eq!(sum, space.vector(&[1.0, 1.0, 0.0]));

    // Scalar multiplication
    let scaled = space.scalar_mul(&v1, 5.0);
    assert_eq!(scaled, space.vector(&[5.0, 0.0, 0.0]));
}

#[test]
fn test_matrix_inverse_laws() {
    for n in 1..=5 {
        let identity = NMatrixMultiplication::zero();
        let mul = NMatrixMultiplication::zero();

        for _ in 0..20 {
            // Create random invertible matrix (for now, diagonal with non-zero entries)
            let elements: Vec<f64> = (0..n).map(|_| rand::random::<f64>().abs() + 1e-1).collect();
            let matrix = NMatrixMultiplication { elements: elements.clone() };

            // Inverse
            let inv = mul.inverse(matrix.clone());
            let product = mul.operate_binary(matrix.clone(), inv.clone());

            // Verify: A * A⁻¹ = I
            for (p, i) in product.elements.iter().zip(identity.elements.iter()) {
                assert!((p - i).abs() < 1e-6, "Matrix inverse product mismatch: {} != {}", p, i);
            }

            // Also verify: A⁻¹ * A = I
            let reverse_product = mul.operate_binary(inv.clone(), matrix.clone());
            for (p, i) in reverse_product.elements.iter().zip(identity.elements.iter()) {
                assert!((p - i).abs() < 1e-6, "Matrix reverse inverse mismatch: {} != {}", p, i);
            }
        }
    }
}

// Matrix group: integer matrix addition
#[test]
fn test_n_matrix_addition_groups() {
    for n in 1..=10 {
        let add = NMatrixAddition::zero();
        let zero = NMatrixAddition::zero();

        for a in 0..n {
            for b in 0..n {
                let mat_a = NMatrixAddition { elements: vec![a as f64; n] };
                let mat_b = NMatrixAddition { elements: vec![b as f64; n] };

                // Closure
                let ab = add.operate_binary(mat_a.clone(), mat_b.clone());
                assert_eq!(ab.elements.len(), n);

                // Commutativity
                assert!(ab.elements.iter().zip(mat_b.elements.iter()).all(|(x, y)| x == y));

                // Associativity
                let ab_c = add.operate_binary(ab.clone(), zero.clone());
                let a_bc = add.operate_binary(mat_a.clone(), add.operate_binary(mat_b.clone(), zero.clone()));
                assert!(ab_c.elements.iter().zip(a_bc.elements.iter()).all(|(x, y)| x == y));

                // Identity
                let identity_result = add.operate_binary(mat_a.clone(), zero.clone());
                assert!(identity_result.elements.iter().zip(mat_a.elements.iter()).all(|(x, y)| x == y));

                // Inverse
                let inv = add.inverse(mat_a.clone());
                let inverse_result = add.operate_binary(mat_a.clone(), inv);
                assert!(inverse_result.elements.iter().zip(zero.elements.iter()).all(|(x, y)| x == y));
            }
        }
    }
}


#[test]
fn test_n_matrix_multiplication_groups() {
    for n in 1..=10 {
        let mul = NMatrixMultiplication::zero();
        let identity = NMatrixMultiplication::zero();

        for a in 1..=n {
            for b in 1..=n {
                let mat_a = NMatrixMultiplication { elements: vec![a as f64; n] };
                let mat_b = NMatrixMultiplication { elements: vec![b as f64; n] };

                // Closure
                let ab = mul.operate_binary(mat_a.clone(), mat_b.clone());
                assert_eq!(ab.elements.len(), n);

                // Associativity
                let ab_c = mul.operate_binary(ab.clone(), identity.clone());
                let a_bc = mul.operate_binary(mat_a.clone(), mul.operate_binary(mat_b.clone(), identity.clone()));
                assert!(ab_c.elements.iter().zip(a_bc.elements.iter()).all(|(x, y)| x == y));

                // Identity
                let identity_result = mul.operate_binary(mat_a.clone(), identity.clone());
                assert!(identity_result.elements.iter().zip(mat_a.elements.iter()).all(|(x, y)| x == y));

                // Inverse
                let inv = mul.inverse(mat_a.clone());
                let inverse_result = mul.operate_binary(mat_a.clone(), inv);
                assert!(inverse_result.elements.iter().zip(identity.elements.iter()).all(|(x, y)| x == y));
            }
        }
    }
}