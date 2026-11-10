file(MAKE_DIRECTORY "${OUTPUT}")
execute_process(COMMAND "${JAVAC}" -cp "${JAR_FILE}" -d "${OUTPUT}" "${SOURCE}" COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${JAVA}" "-Djava.library.path=${LIB_DIR}" -cp "${JAR_FILE}:${OUTPUT}" Smoke COMMAND_ERROR_IS_FATAL ANY)
