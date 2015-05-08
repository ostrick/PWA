FILE(REMOVE_RECURSE
  "CMakeFiles/Build"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang)
  INCLUDE(CMakeFiles/Build.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
