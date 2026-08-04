@echo off
if defined RASPA3_EXECUTABLE (
  "%RASPA3_EXECUTABLE%" %*
) else (
  raspa3.exe %*
)
