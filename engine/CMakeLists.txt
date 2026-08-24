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
    src/matrix/tms_local_matrix.f90
    src/geometry/tms_geometry.f90
    src/materials/tms_dynamic_modulus.f90
    src/materials/tms_material_frequency.f90
    src/materials/tms_material.f90
    src/solver/tms_inertia.f90
    src/solver/tms_torsional_stiffness.f90
    src/solver/tms_dynamic_torsional_stiffness.f90
    src/solver/tms_frequency_solver.f90
    src/system/tms_torsional_node.f90
    src/system/tms_torsional_element.f90
    src/system/tms_torsional_system.f90
    src/solver/tms_torsional_system.f90
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

if(BUILD_TESTING)
  # Tüm Fortran test hedeflerine aynı standart ve bağlantı kurallarını uygular.
  function(tms26_add_fortran_test test_name source_file)
    set(test_target "tms26_test_${test_name}")

    add_executable("${test_target}" "${source_file}")
    target_link_libraries("${test_target}" PRIVATE tms26_engine)

    set_target_properties(
      "${test_target}"
      PROPERTIES
        Fortran_STANDARD 2018
        Fortran_STANDARD_REQUIRED YES
        Fortran_EXTENSIONS NO
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
endif()
