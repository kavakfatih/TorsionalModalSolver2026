cmake_minimum_required(VERSION 3.25)

project(
  TorsionalModalSolver2026
  VERSION 0.2.2
  DESCRIPTION "Elastomer esaslı burulma titreşimi sistemleri hesap motoru"
  LANGUAGES Fortran
)

# IDE ve kod inceleme araçlarının derleme komutlarına erişmesini sağlar.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# BUILD_TESTING seçeneğini tanımlar ve test altyapısını etkinleştirir.
include(CTest)

add_subdirectory(engine)
