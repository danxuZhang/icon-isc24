## Formatting
clang-format is used to enforce a consistent code style. A .clang-format file is provided in the repository.
To format all files in place you can use:

`clang-format -i -style=file core/**/*.cpp core/**/*.hpp implementations/**/*.cpp implementations/**/*.hpp io/**/*.cpp io/**/*.hpp test/*.cc main.cpp`.

Because development concentrates mainly on main.cpp and the implementations, this shortened command can also be used: 

`clang-format -i -style=file implementations/**/*.cpp implementations/**/*.hpp test/*.cc main.cpp`.

We have to specify the folders because e.g. `seq` will build their files in a subfolder as well and we don't want to format their files.

It is recommended to create a precommit file for Git so the changes are formatted automatically. For this we have to create .git/hook/pre-commit.
We also need to make sure that clang-format is installed.
A simple example would be:
`#!/bin/bash
echo Running clang-format on added files.
for FILE in $(git diff --cached --name-only | grep -E '.*\.(c|cpp|h|hpp)\b')
do
  clang-format -i $FILE
done
`
Make sure to mark the file as executable with e.g. `chmod +x .git/hook/pre-commit`.

