# TMS26 hesap motorunun kaynaklarını tek bir statik kitaplıkta toplar.
add_library(tms26_engine STATIC)

target_sources(
  tms26_engine
  PRIVATE
    # Fortran modülleri sağlayıcıdan tüketiciye doğru listelenir. Gerçek
    # derleme bağımlılıkları CMake Fortran tarayıcısı tarafından oluşturulur.
    src/core/tms_kinds.f90
    src/core/tms_constants.f90
    src/core/tms_units.f90
    src/numerics/tms_scalar_minimizer.f90
    src/numerics/tms_weighted_piecewise_integrals.f90
    src/numerics/tms_sample_statistics.f90
    src/numerics/tms_deterministic_rng.f90
    src/numerics/tms_bootstrap.f90
    src/identification/tms_tts_types.f90
    src/identification/tms_tts_uncertainty_types.f90
    src/identification/tms_tts_shift_law_types.f90
    src/identification/tms_tts_shift_law_fit.f90
    src/identification/tms_tts_shift_law_validation.f90
    src/identification/tms_tts_repeatability_types.f90
    src/identification/tms_tts_repeatability_analysis.f90
    src/identification/tms_tts_pair_shift.f90
    src/identification/tms_tts_uncertainty_validation.f90
    src/identification/tms_tts_weighted_pair_shift.f90
    src/identification/tms_tts_uncertainty_sensitivity.f90
    src/identification/tms_tts_shift_chain.f90
    src/identification/tms_tts_master_curve.f90
    src/identification/tms_tts_diagnostics.f90
    src/system/tms_dof_types.f90
    src/matrix/tms_matrix_types.f90
    src/matrix/tms_local_matrix.f90
    src/geometry/tms_geometry.f90
    src/materials/tms_dynamic_modulus.f90
    src/materials/tms_material_frequency.f90
    src/materials/tms_dynamic_material_metadata.f90
    src/materials/tms_temperature_shift_types.f90
    src/materials/tms_dynamic_modulus_provider.f90
    src/materials/tms_temperature_shift_provider.f90
    src/materials/tms_tabulated_dynamic_modulus_provider.f90
    src/materials/tms_tabulated_temperature_shift.f90
    src/materials/tms_wlf_temperature_shift.f90
    src/materials/tms_arrhenius_temperature_shift.f90
    src/materials/tms_thermorheological_dynamic_modulus_provider.f90
    src/identification/tms_tts_identification.f90
    src/identification/tms_tts_runtime_export.f90
    src/identification/tms_tts_parametric_runtime_export.f90
    src/materials/tms_material.f90
    src/solver/tms_inertia.f90
    src/solver/tms_torsional_stiffness.f90
    src/solver/tms_dynamic_torsional_stiffness.f90
    src/solver/tms_frequency_solver.f90
    src/system/tms_dynamic_torsional_property_binding.f90
    src/system/tms_torsional_node.f90
    src/system/tms_torsional_element.f90
    src/system/tms_torsional_system.f90
    src/matrix/tms_dof_map.f90
    src/constraints/tms_constraint_types.f90
    src/constraints/tms_constraint_manager.f90
    src/matrix/tms_stiffness_matrix.f90
    src/matrix/tms_loss_stiffness_matrix.f90
    src/matrix/tms_damping_matrix.f90
    src/matrix/tms_mass_matrix.f90
    src/matrix/tms_matrix_reduction.f90
    src/matrix/tms_matrix_assembly.f90
    src/constraints/tms_reduced_system.f90
    src/harmonic/tms_reduced_dynamic_system.f90
    src/solver/tms_torsional_system.f90
    src/eigen/tms_generalized_eigen_problem.f90
    src/eigen/tms_eigen_solution.f90
    src/eigen/tms_lapack_dsygv_backend.f90
    src/eigen/tms_generalized_eigen_solver.f90
    src/modal/tms_modal_validation.f90
    src/modal/tms_modal_result.f90
    src/modal/tms_modal_analysis.f90
    src/harmonic/tms_complex_linear_problem.f90
    src/harmonic/tms_complex_linear_solution.f90
    src/harmonic/tms_lapack_zsysvx_backend.f90
    src/harmonic/tms_complex_linear_solver.f90
    src/harmonic/tms_dynamic_stiffness.f90
    src/harmonic/tms_harmonic_excitation.f90
    src/harmonic/tms_harmonic_response.f90
    src/harmonic/tms_frf.f90
    src/harmonic/tms_material_state_trace.f90
    src/harmonic/tms_material_aware_harmonic_response.f90
    src/harmonic/tms_harmonic_analysis.f90
    src/harmonic/tms_material_aware_harmonic_analysis.f90
)

set_target_properties(
  tms26_engine
  PROPERTIES
    Fortran_STANDARD 2018
    Fortran_STANDARD_REQUIRED YES
    Fortran_EXTENSIONS NO
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
)

# Bağlanan hedeflerin üretilen Fortran modüllerini bulabilmesini sağlar.
target_include_directories(
  tms26_engine
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>"
)

# DSYGV backend'inin platform bağımsız LAPACK bağlantısını sağlar.
target_link_libraries(tms26_engine PRIVATE LAPACK::LAPACK)

if(BUILD_TESTING)
  # Tüm Fortran test hedeflerine aynı standart ve bağlantı kurallarını uygular.
  function(tms26_add_fortran_test test_name source_file)
    set(test_target "tms26_test_${test_name}")

    add_executable("${test_target}" "${source_file}" ${ARGN})
    target_link_libraries("${test_target}" PRIVATE tms26_engine)

    set_target_properties(
      "${test_target}"
      PROPERTIES
        Fortran_STANDARD 2018
        Fortran_STANDARD_REQUIRED YES
        Fortran_EXTENSIONS NO
        Fortran_MODULE_DIRECTORY
          "${CMAKE_CURRENT_BINARY_DIR}/test-modules/${test_target}"
    )

    add_test(NAME "tms26.${test_name}" COMMAND "${test_target}")
  endfunction()

  tms26_add_fortran_test(kinds tests/test_kinds.f90)
  tms26_add_fortran_test(constants tests/test_constants.f90)
  tms26_add_fortran_test(units tests/test_units.f90)
  tms26_add_fortran_test(geometry tests/test_geometry.f90)
  tms26_add_fortran_test(material tests/test_material.f90)
  tms26_add_fortran_test(dynamic_modulus tests/test_dynamic_modulus.f90)
  tms26_add_fortran_test(torsional_node tests/test_torsional_node.f90)
  tms26_add_fortran_test(torsional_element tests/test_torsional_element.f90)
  tms26_add_fortran_test(
    local_stiffness_matrix
    tests/test_local_stiffness_matrix.f90
  )
  add_test(
    NAME "tms26.local_stiffness_matrix.rejects_negative_stiffness"
    COMMAND
      "$<TARGET_FILE:tms26_test_local_stiffness_matrix>"
      "negative_stiffness"
  )
  set_tests_properties(
    "tms26.local_stiffness_matrix.rejects_negative_stiffness"
    PROPERTIES WILL_FAIL TRUE
  )
  tms26_add_fortran_test(
    generalized_torsional_system
    tests/test_generalized_torsional_system.f90
  )
  tms26_add_fortran_test(matrix_assembly tests/test_matrix_assembly.f90)
  tms26_add_fortran_test(
    constraint_foundation
    tests/test_constraint_foundation.f90
  )

  foreach(
    validation_case
    IN ITEMS
      missing_node
      duplicate_physical_dof
      duplicate_constraint_id
      invalid_dof_type
      invalid_constraint_type
      nonfinite_value
      fixed_nonzero_value
      legacy_constraint_mismatch
      foreign_full_dof_map
      invalid_recovery_size
      invalid_recovery_value
  )
    set(
      validation_test
      "tms26.constraint_foundation.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_constraint_foundation>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    torsional_validation
    tests/test_torsional_validation.f90
  )

  foreach(
    validation_case
    IN ITEMS
      nan_tolerance
      positive_infinity_tolerance
      negative_infinity_tolerance
  )
    set(
      validation_test
      "tms26.torsional_validation.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_torsional_validation>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    numerical_hardening
    tests/test_numerical_hardening.f90
  )

  foreach(
    nonfinite_kind
    IN ITEMS
      nan
      positive_infinity
      negative_infinity
  )
    foreach(
      physical_quantity
      IN ITEMS
        polar_inertia
        density
        radius
        length
        stiffness
        modulus
        frequency
        temperature
    )
      set(validation_case "${nonfinite_kind}_${physical_quantity}")
      set(
        validation_test
        "tms26.numerical_hardening.rejects_${validation_case}"
      )
      add_test(
        NAME "${validation_test}"
        COMMAND
          "$<TARGET_FILE:tms26_test_numerical_hardening>"
          "${validation_case}"
      )
      set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
    endforeach()
  endforeach()

  foreach(
    validation_case
    IN ITEMS
      hub_rubber_interface_mismatch
      rubber_ring_interface_mismatch
      zero_length
      reversed_annular_geometry
      output_overflow_geometry
      output_overflow_inertia
      output_overflow_stiffness
      output_overflow_loss_factor
  )
    set(
      validation_test
      "tms26.numerical_hardening.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_numerical_hardening>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  foreach(
    validation_case
    IN ITEMS
      negative_matrix_dimension
      unknown_node_lookup
      out_of_range_equation
      nonpositive_inertia
      incompatible_dof_map
  )
    set(validation_test "tms26.matrix_assembly.rejects_${validation_case}")
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_matrix_assembly>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  foreach(
    validation_case
    IN ITEMS
      nonpositive_id
      nonpositive_inertia
      nonfinite_angle
  )
    set(validation_test "tms26.torsional_node.rejects_${validation_case}")
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_torsional_node>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  foreach(
    validation_case
    IN ITEMS
      nonpositive_id
      nonpositive_node_id
      self_connection
      nonpositive_stiffness
      negative_damping
      nonfinite_stiffness
  )
    set(validation_test "tms26.torsional_element.rejects_${validation_case}")
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_torsional_element>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  foreach(
    validation_case
    IN ITEMS
      empty_system
      duplicate_node_id
      missing_endpoint
      duplicate_element_id
  )
    set(
      validation_test
      "tms26.generalized_torsional_system.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_generalized_torsional_system>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(inertia tests/test_inertia.f90)
  tms26_add_fortran_test(hub_inertia tests/test_hub_inertia.f90)
  tms26_add_fortran_test(torsional_stiffness tests/test_torsional_stiffness.f90)
  add_test(
    NAME "tms26.torsional_stiffness.rejects_nonpositive_modulus"
    COMMAND
      "$<TARGET_FILE:tms26_test_torsional_stiffness>"
      "nonpositive_modulus"
  )
  set_tests_properties(
    "tms26.torsional_stiffness.rejects_nonpositive_modulus"
    PROPERTIES WILL_FAIL TRUE
  )
  tms26_add_fortran_test(
    dynamic_torsional_stiffness
    tests/test_dynamic_torsional_stiffness.f90
  )

  # Aynı yürütülebilir dosyayı geçersiz fiziksel girdilerle çalıştırır.
  # Üretim yordamının error stop ile reddetmesi bu testlerin başarı koşuludur.
  foreach(
    validation_case
    IN ITEMS
      negative_radius
      zero_inner_radius
      unordered_radii
      equal_radii
      zero_length
      negative_length
      nonpositive_storage_modulus
      negative_loss_modulus
  )
    set(
      validation_test
      "tms26.dynamic_torsional_stiffness.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_dynamic_torsional_stiffness>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(frequency_solver tests/test_frequency_solver.f90)
  tms26_add_fortran_test(torsional_system tests/test_torsional_system.f90)

  tms26_add_fortran_test(
    generalized_eigen_solver
    tests/test_generalized_eigen_solver.f90
  )

  # Generalized eigenproblem structural contract'ı ile DSYGV M-SPD tanısını
  # ayrı süreçlerde sınar. Üretim API'sinin error stop vermesi başarıdır.
  foreach(
    validation_case
    IN ITEMS
      nonsquare_stiffness
      nonsquare_mass
      dimension_mismatch
      nonsymmetric_stiffness
      nonsymmetric_mass
      nan_stiffness
      positive_infinity_stiffness
      negative_infinity_stiffness
      nan_mass
      positive_infinity_mass
      negative_infinity_mass
      zero_order
      excessive_mode_count
  )
    set(
      validation_test
      "tms26.generalized_eigen_solver.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_generalized_eigen_solver>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  # DSYGV INFO>N yolunun yalnız başarısız çıkmasını değil, raw INFO yerine
  # anlamlı TMS26 M-SPD tanısını üretmesini de platform bağımsız doğrular.
  foreach(
    validation_case
    IN ITEMS
      zero_mass
      negative_mass
      indefinite_mass
  )
    set(
      validation_test
      "tms26.generalized_eigen_solver.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "${CMAKE_COMMAND}"
        "-DTEST_EXECUTABLE=$<TARGET_FILE:tms26_test_generalized_eigen_solver>"
        "-DTEST_SELECTOR=${validation_case}"
        "-DEXPECTED_REGEX=M positive definite"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/assert_expected_failure.cmake"
    )
  endforeach()

  tms26_add_fortran_test(modal_analysis tests/test_modal_analysis.f90)

  # Seçenek hata yollarını ayrı CTest süreçlerinde doğrular.
  foreach(
    validation_case
    IN ITEMS
      negative_requested_mode_count
      excessive_requested_mode_count
      nonfinite_rigid_multiplier
  )
    set(validation_test "tms26.modal_analysis.rejects_${validation_case}")
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_modal_analysis>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  # Fully constrained ve significant negative eigenvalue vakalarının doğru
  # pre/post-processing tanı yolunda, LAPACK misuse olmadan reddini kilitler.
  foreach(validation_case IN ITEMS fully_constrained significant_negative_eigenvalue)
    if(validation_case STREQUAL "fully_constrained")
      set(expected_regex "aktif DOF yoktur")
    else()
      set(expected_regex "negatif")
    endif()
    set(validation_test "tms26.modal_analysis.rejects_${validation_case}")
    add_test(
      NAME "${validation_test}"
      COMMAND
        "${CMAKE_COMMAND}"
        "-DTEST_EXECUTABLE=$<TARGET_FILE:tms26_test_modal_analysis>"
        "-DTEST_SELECTOR=${validation_case}"
        "-DEXPECTED_REGEX=${expected_regex}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/tests/assert_expected_failure.cmake"
    )
  endforeach()

  # Sistem katmanının doğrudan modal girdileri ile builder yoğunluk ve
  # geometri önkoşullarını error stop üzerinden reddettiğini doğrular.
  foreach(
    validation_case
    IN ITEMS
      nonpositive_hub_inertia
      nonpositive_ring_inertia
      nonpositive_storage_stiffness
      nonpositive_hub_density
      nonpositive_ring_density
      invalid_hub_geometry
      invalid_ring_geometry
  )
    set(validation_test "tms26.torsional_system.rejects_${validation_case}")
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_torsional_system>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    dynamic_matrix_assembly
    tests/test_dynamic_matrix_assembly.f90
  )

  foreach(
    validation_case
    IN ITEMS
      negative_loss_stiffness
      nan_loss_stiffness
      positive_infinity_loss_stiffness
      negative_infinity_loss_stiffness
      negative_damping
      nan_damping
      positive_infinity_damping
      negative_infinity_damping
  )
    set(
      validation_test
      "tms26.dynamic_matrix_assembly.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_dynamic_matrix_assembly>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    complex_linear_solver
    tests/test_complex_linear_solver.f90
  )

  foreach(
    validation_case
    IN ITEMS
      zero_order
      zero_rhs
      nonsquare_coefficient
      rhs_dimension_mismatch
      nonsymmetric_coefficient
      nan_coefficient_real
      infinite_coefficient_imaginary
      nan_rhs_imaginary
      infinite_rhs_real
      singular_response_access
  )
    set(
      validation_test
      "tms26.complex_linear_solver.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_complex_linear_solver>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    harmonic_analysis
    tests/test_harmonic_analysis.f90
  )

  foreach(
    validation_case
    IN ITEMS
      empty_frequency
      zero_frequency
      negative_frequency
      nan_frequency
      positive_infinity_frequency
      negative_infinity_frequency
      unordered_frequency
      duplicate_frequency
      empty_excitation
      unknown_node
      unsupported_dof
      constrained_target
      nan_torque_real
      nan_torque_imaginary
      infinite_torque_real
      infinite_torque_imaginary
      overflow_torque_sum
      fully_constrained
      zero_frf_torque
  )
    set(
      validation_test
      "tms26.harmonic_analysis.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_harmonic_analysis>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    tabulated_dynamic_modulus_provider
    tests/test_tabulated_dynamic_modulus_provider.f90
  )

  foreach(
    validation_case
    IN ITEMS
      less_than_two_points
      zero_frequency
      negative_frequency
      nan_frequency
      infinite_frequency
      duplicate_frequency
      numerical_duplicate_frequency
      unordered_frequency
      zero_storage_modulus
      negative_storage_modulus
      nan_storage_modulus
      infinite_storage_modulus
      negative_loss_modulus
      nan_loss_modulus
      infinite_loss_modulus
      invalid_dataset_temperature
      nan_dataset_temperature
      infinite_dataset_temperature
      mixed_point_temperature
      invalid_policy
      below_linear_domain
      above_linear_domain
      below_log_domain
      above_log_domain
      temperature_mismatch
      invalid_strain_amplitude
      invalid_prestrain
      magic_unavailable_strain
      magic_unavailable_prestrain
      invalid_mode
  )
    set(
      validation_test
      "tms26.tabulated_dynamic_modulus_provider.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_tabulated_dynamic_modulus_provider>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    dynamic_torsional_property_binding
    tests/test_dynamic_torsional_property_binding.f90
  )

  foreach(
    validation_case
    IN ITEMS
      nonpositive_element_id
      negative_radius
      unordered_radii
      zero_length
      uninitialized_provider
      tensile_mode
  )
    set(
      validation_test
      "tms26.dynamic_torsional_property_binding.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_dynamic_torsional_property_binding>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    material_aware_harmonic_analysis
    tests/test_material_aware_harmonic_analysis.f90
  )

  foreach(
    validation_case
    IN ITEMS
      duplicate_binding
      unknown_element
      uncovered_sweep
      temperature_mismatch
      empty_bindings
  )
    set(
      validation_test
      "tms26.material_aware_harmonic_analysis.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_material_aware_harmonic_analysis>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    tabulated_temperature_shift
    tests/test_tabulated_temperature_shift.f90
  )

  foreach(
    validation_case
    IN ITEMS
      less_than_two_points
      size_mismatch
      zero_temperature
      negative_temperature
      nan_temperature
      infinite_temperature
      duplicate_temperature
      numerical_duplicate_temperature
      unordered_temperature
      nan_shift
      infinite_shift
      invalid_reference_temperature
      missing_reference_temperature
      nonzero_reference_shift
      below_domain
      above_domain
      inconsistent_shift_factor
      invalid_shift_model_kind
      invalid_temperature_bracket
      nonfinite_temperature_bracket
      inconsistent_provider_evaluation
  )
    set(
      validation_test
      "tms26.tabulated_temperature_shift.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_tabulated_temperature_shift>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    wlf_temperature_shift
    tests/test_wlf_temperature_shift.f90
  )

  foreach(
    validation_case
    IN ITEMS
      zero_c1
      negative_c1
      nan_c1
      zero_c2
      negative_c2
      infinite_c2
      invalid_reference_temperature
      invalid_minimum_temperature
      reversed_domain
      reference_outside_domain
      pole_in_domain
      below_domain
      above_domain
      invalid_query_temperature
  )
    set(
      validation_test
      "tms26.wlf_temperature_shift.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_wlf_temperature_shift>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    arrhenius_temperature_shift
    tests/test_arrhenius_temperature_shift.f90
  )

  foreach(
    validation_case
    IN ITEMS
      zero_activation_energy
      negative_activation_energy
      nan_activation_energy
      infinite_activation_energy
      invalid_reference_temperature
      invalid_minimum_temperature
      reversed_domain
      reference_outside_domain
      below_domain
      above_domain
      zero_query_temperature
      nan_query_temperature
  )
    set(
      validation_test
      "tms26.arrhenius_temperature_shift.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_arrhenius_temperature_shift>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    thermorheological_dynamic_modulus_provider
    tests/test_thermorheological_dynamic_modulus_provider.f90
  )

  foreach(
    validation_case
    IN ITEMS
      reference_temperature_mismatch
      temperature_below_domain
      temperature_above_domain
      reduced_frequency_below_domain
      reduced_frequency_above_domain
      zero_physical_frequency
      nan_physical_frequency
      infinite_physical_frequency
      logspace_overflow_outside_domain
      logspace_underflow_outside_domain
  )
    set(
      validation_test
      "tms26.thermorheological_dynamic_modulus_provider.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_thermorheological_dynamic_modulus_provider>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  tms26_add_fortran_test(
    thermorheological_harmonic_analysis
    tests/test_thermorheological_harmonic_analysis.f90
  )

  foreach(
    validation_case
    IN ITEMS
      uncovered_reduced_sweep
      temperature_outside_domain
      reference_trace_metadata_mismatch
  )
    set(
      validation_test
      "tms26.thermorheological_harmonic_analysis.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_thermorheological_harmonic_analysis>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  set(tts_test_support tests/tms_tts_test_support.f90)
  tms26_add_fortran_test(
    tts_data_model
    tests/test_tts_data_model.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_scalar_minimizer
    tests/test_tts_scalar_minimizer.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_pair_shift
    tests/test_tts_pair_shift.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_shift_chain
    tests/test_tts_shift_chain.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_master_curve
    tests/test_tts_master_curve.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_diagnostics
    tests/test_tts_diagnostics.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_identification
    tests/test_tts_identification.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_generalized_maxwell
    tests/test_tts_generalized_maxwell.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_runtime_roundtrip
    tests/test_tts_runtime_roundtrip.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_arrhenius_fit
    tests/test_tts_arrhenius_fit.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_wlf_fit
    tests/test_tts_wlf_fit.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_shift_law_validation
    tests/test_tts_shift_law_validation.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_parametric_runtime_roundtrip
    tests/test_tts_parametric_runtime_roundtrip.f90
    ${tts_test_support}
  )
  foreach(
    validation_case
    IN ITEMS
      arrhenius_outside_domain
      wlf_outside_domain
  )
    set(
      validation_test
      "tms26.tts_parametric_runtime_roundtrip.rejects_${validation_case}"
    )
    add_test(
      NAME "${validation_test}"
      COMMAND
        "$<TARGET_FILE:tms26_test_tts_parametric_runtime_roundtrip>"
        "${validation_case}"
    )
    set_tests_properties("${validation_test}" PROPERTIES WILL_FAIL TRUE)
  endforeach()

  set(
    tts_repeatability_test_support
    tests/tms_tts_repeatability_test_support.f90
  )
  set(
    tts_uncertainty_test_support
    tests/tms_tts_uncertainty_test_support.f90
  )
  tms26_add_fortran_test(
    tts_sample_statistics
    tests/test_tts_sample_statistics.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_bootstrap
    tests/test_tts_bootstrap.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_repeatability
    tests/test_tts_repeatability.f90
    ${tts_test_support}
    ${tts_repeatability_test_support}
  )
  tms26_add_fortran_test(
    tts_parametric_repeatability
    tests/test_tts_parametric_repeatability.f90
    ${tts_test_support}
    ${tts_repeatability_test_support}
  )
  tms26_add_fortran_test(
    tts_repeatability_bootstrap
    tests/test_tts_repeatability_bootstrap.f90
    ${tts_test_support}
    ${tts_repeatability_test_support}
  )
  tms26_add_fortran_test(
    tts_uncertainty_propagation
    tests/test_tts_uncertainty_propagation.f90
    ${tts_test_support}
    ${tts_uncertainty_test_support}
  )
  tms26_add_fortran_test(
    tts_weighted_integrals
    tests/test_tts_weighted_integrals.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_huber_integrals
    tests/test_tts_huber_integrals.f90
    ${tts_test_support}
  )
  tms26_add_fortran_test(
    tts_weighted_pair_shift
    tests/test_tts_weighted_pair_shift.f90
    ${tts_test_support}
    ${tts_uncertainty_test_support}
  )
  tms26_add_fortran_test(
    tts_robust_shift_sensitivity
    tests/test_tts_robust_shift_sensitivity.f90
    ${tts_test_support}
    ${tts_uncertainty_test_support}
  )
endif()
