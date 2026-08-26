# Bir üretim executable'ının beklenen selector ile başarısız olduğunu ve
# tanısal çıktısının gerekli ifadeyi içerdiğini doğrular. Bu sarmalayıcı,
# CTest WILL_FAIL özelliğinin yalnız exit code denetimiyle yetinmediği hata
# sözleşmelerinde platform bağımsız olarak kullanılır.

foreach(required_variable IN ITEMS TEST_EXECUTABLE TEST_SELECTOR EXPECTED_REGEX)
  if(NOT DEFINED ${required_variable})
    message(FATAL_ERROR "Eksik test sarmalayıcı girdisi: ${required_variable}")
  endif()
endforeach()

execute_process(
  COMMAND "${TEST_EXECUTABLE}" "${TEST_SELECTOR}"
  RESULT_VARIABLE test_result
  OUTPUT_VARIABLE test_stdout
  ERROR_VARIABLE test_stderr
)

if("${test_result}" STREQUAL "0")
  message(FATAL_ERROR "Geçersiz vaka beklenmedik biçimde başarılı oldu.")
endif()

string(CONCAT test_output "${test_stdout}" "\n" "${test_stderr}")
if(NOT "${test_output}" MATCHES "${EXPECTED_REGEX}")
  message(FATAL_ERROR
    "Beklenen tanı bulunamadı: ${EXPECTED_REGEX}\nGerçek çıktı:\n${test_output}"
  )
endif()
