#[test]
fn test_ndim_number_operations() {
    let add = NDimAddition;

    let a = NDimNumber::new(vec![1.0, 2.0, 3.0]);
    let b = NDimNumber::new(vec![4.0, -5.0, 6.0]);

    // Closure
    let result = add.operate_binary(a.clone(), b.clone());
    assert_eq!(result.components, vec![5.0, -3.0, 9.0]);

    // Identity
    let identity = NDimNumber::zero(3);
    assert_eq!(add.operate_binary(a.clone(), identity.clone()), a);

    // Inverse
    let inverse = add.inverse(a.clone());
    assert_eq!(add.operate_binary(a.clone(), inverse), identity);
}

#[test]
fn test_euclidean_space_vector_operations() {
    let space = EuclideanSpace::new(3);

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
fn test_nmatrix_addition_operations() {
    let add = NMatrixAddition::zero();
    let zero = NMatrixAddition::zero();

    let a = NMatrixAddition { elements: vec![1.0, 2.0] };
    let b = NMatrixAddition { elements: vec![3.0, 4.0] };

    // Closure
    let result = add.operate_binary(a.clone(), b.clone());
    assert_eq!(result.elements, vec![4.0, 6.0]);

    // Identity
    let identity_result = add.operate_binary(a.clone(), zero.clone());
    assert_eq!(identity_result.elements, a.elements);

    // Inverse
    let inverse = add.inverse(a.clone());
    let inverse_result = add.operate_binary(a.clone(), inverse);
    assert_eq!(inverse_result.elements, zero.elements);
}

#[test]
fn test_nmatrix_multiplication_operations() {
    let mul = NMatrixMultiplication::zero();
    let identity = NMatrixMultiplication::zero();

    let a = NMatrixMultiplication { elements: vec![2.0, 3.0] };
    let b = NMatrixMultiplication { elements: vec![4.0, 5.0] };

    // Closure
    let result = mul.operate_binary(a.clone(), b.clone());
    assert_eq!(result.elements, vec![8.0, 15.0]);

    // Identity
    let identity_result = mul.operate_binary(a.clone(), identity.clone());
    assert_eq!(identity_result.elements, a.elements);

    // Inverse
    let inverse = mul.inverse(a.clone());
    let inverse_result = mul.operate_binary(a.clone(), inverse);
    assert_eq!(inverse_result.elements, identity.elements);
}

#[test]
fn test_quaternion_norm_and_inverse() {
    let mul = QuaternionMultiplication;

    let a = Quaternion { re: 1.0, i: 2.0, j: 3.0, k: 4.0 };

    // Norm
    let norm_sq = a.re * a.re + a.i * a.i + a.j * a.j + a.k * a.k;
    assert_eq!(norm_sq, 30.0);

    // Inverse
    let inv = mul.inverse(a);
    let product = mul.operate_binary(a, inv);
    let identity = mul.identity();
    assert_eq!(product, identity);
}

#[test]
#[should_panic(expected = "Dimensions must match")]
fn test_ndim_addition_mismatched_dimensions() {
    let add = NDimAddition;

    let a = NDimNumber::new(vec![1.0, 2.0]);
    let b = NDimNumber::new(vec![3.0]);

    add.operate_binary(a, b); // Should panic
}

#[test]
fn test_empty_input_operations() {
    let add = NDimAddition;
    let mul = NMatrixMultiplication::zero();

    // Empty addition
    let result_add = add.operate(&[]);
    assert_eq!(result_add.components, vec![]);

    // Empty multiplication
    let result_mul = mul.operate(&[]);
    assert_eq!(result_mul.elements, vec![1.0; 2]); // Identity for multiplication
}

#[test]
fn test_trait_implementation() {
    let add = NDimAddition;
    let a = NDimNumber::new(vec![1.0, 2.0, 3.0]);
    let b = NDimNumber::new(vec![4.0, 5.0, 6.0]);

    // Test Algebra trait
    let result = add.operate_binary(a.clone(), b.clone());
    assert_eq!(result.components, vec![5.0, 7.0, 9.0]);

    // Test Group trait
    let identity = add.identity();
    assert_eq!(identity.components, vec![0.0, 0.0, 0.0]);

    let inverse = add.inverse(a.clone());
    assert_eq!(inverse.components, vec![-1.0, -2.0, -3.0]);
}