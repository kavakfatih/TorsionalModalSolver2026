# TMS26 hesap motorunun kaynaklarını tek bir statik kitaplıkta toplar.
add_library(tms26_engine STATIC)

target_sources(
  tms26_engine
  PRIVATE
    src/core/tms_kinds.f90
    src/core/tms_constants.f90
    src/core/tms_units.f90
    src/geometry/tms_geometry.f90
    src/materials/tms_dynamic_modulus.f90
    src/materials/tms_material_frequency.f90
    src/materials/tms_material.f90
    src/solver/tms_inertia.f90
    src/solver/tms_torsional_stiffness.f90
    src/solver/tms_dynamic_torsional_stiffness.f90
    src/solver/tms_frequency_solver.f90
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
  tms26_add_fortran_test(inertia tests/test_inertia.f90)
  tms26_add_fortran_test(torsional_stiffness tests/test_torsional_stiffness.f90)
  tms26_add_fortran_test(
    dynamic_torsional_stiffness
    tests/test_dynamic_torsional_stiffness.f90
  )
  tms26_add_fortran_test(frequency_solver tests/test_frequency_solver.f90)
endif()
