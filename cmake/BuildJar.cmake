file(GLOB sources "${JAVA_DIR}/*.java")
file(MAKE_DIRECTORY "${JAVA_DIR}/classes")
execute_process(COMMAND "${JAVAC}" --release 11 -d "${JAVA_DIR}/classes" ${sources} COMMAND_ERROR_IS_FATAL ANY)
execute_process(COMMAND "${JAR}" --create --file "${OUTPUT}" -C "${JAVA_DIR}/classes" . COMMAND_ERROR_IS_FATAL ANY)
