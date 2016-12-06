project(RegulationIntegrationTests)

configure_file(run_integration_tests.sh.in "${CMAKE_CURRENT_BINARY_DIR}/run_integration_tests.sh")
configure_file(integration_tests_variables.sh.in "${CMAKE_CURRENT_BINARY_DIR}/integration_tests_variables.sh")

add_custom_target(regulation_integration_test
	COMMAND bash "${CMAKE_CURRENT_BINARY_DIR}/run_integration_tests.sh"
)
