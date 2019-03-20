#include "../IO/H5IODataTree.h"

int main(int, char **)
{
    IO::H5IODataTree tree;
    int test_ints[] = {1, 2, 3, 9, 8, 7};
    double test_doubles[] = {1.0, M_PI, M_PI*M_PI};
    char *test_strings[4];
    for(unsigned i = 0; i < 4; ++i) {
        test_strings[i] = new char[5];
        strcpy(test_strings[i], "str1");
    }
    tree.add_c_array(
        "test.int",
        test_ints,
        "int",
        6
    );
    tree.add_c_array(
        "test.double",
        test_doubles,
        "double",
        3
    );
    tree.add_c_array(
        "test.string",
        test_strings,
        "str",
        4
    );
    for(unsigned i = 0; i < 4; ++i)
        delete[] test_strings[i];
}
