git_root()
{
    git rev-parse --show-toplevel
}

find_source_files()
{
    cd $(git_root)
    pwd
    for d in ./include ./src ./lib ./lib-py; do
        find $d -type f -name "*.hpp" -o -name "*.cpp" -o -name "*.cxx"
    done
}

for f in $(find find_source_files); do
    clang-format-10 -i $f
done
