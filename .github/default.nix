with import <unstable> {};

mkShell {
    buildInputs = [ gcc gdb 
        cppcheck clangAnalyzer clang-tools
    ];
}
