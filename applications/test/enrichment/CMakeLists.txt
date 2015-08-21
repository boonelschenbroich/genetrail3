project(EnrichmentIntegrationTests)

configure_file(run_integration_tests.sh.in "${CMAKE_CURRENT_BINARY_DIR}/run_integration_tests.sh")
configure_file(integration_tests_variables.sh.in "${CMAKE_CURRENT_BINARY_DIR}/integration_tests_variables.sh")

configure_file(input/categories.txt.in "${CMAKE_CURRENT_BINARY_DIR}/input/categories.txt")

add_custom_target(integration_test
	COMMAND bash "${CMAKE_CURRENT_BINARY_DIR}/run_integration_tests.sh"
)
