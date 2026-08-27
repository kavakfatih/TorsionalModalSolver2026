cmake_minimum_required(VERSION 3.25)

project(
  TorsionalModalSolver2026
  VERSION 0.6.0
  DESCRIPTION "Elastomer esaslı burulma titreşimi sistemleri hesap motoru"
  LANGUAGES Fortran
)

# IDE ve kod inceleme araçlarının derleme komutlarına erişmesini sağlar.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# LAPACK ile kullanılan Fortran INTEGER arayüzünü LP64 olarak sabitler.
set(BLA_SIZEOF_INTEGER 4)
find_package(LAPACK REQUIRED)

# BUILD_TESTING seçeneğini tanımlar ve test altyapısını etkinleştirir.
include(CTest)

add_subdirectory(engine)
