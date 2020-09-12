git_root()
{
    git rev-parse --show-toplevel
}

find_source_files()
{
    cd $(git_root)
    for d in ./include ./src ./lib ./lib-py ./tests; do
        find $d -type f -name "*.hpp" -o -name "*.cpp" -o -name "*.cxx"
    done
}

find_py_source_files(){
    cd $(git_root)
    for d in ./lib-py ./tests/py/; do
        find $d -type f -name "*.py" -o -name "*.in"
    done
}

for f in $(find_source_files); do
    echo "Formatting $f ..."
    clang-format-10 -i $f
done

for p in $(find_py_source_files); do
    echo "Formatting $p ..."
    yapf -i --style='{based_on_style:pep8}' $p
done

