#include "../IO/H5IODataTree.h"

int main(int, char **)
{
    IO::H5IODataTree tree;
    int test_ints[] = {1, 2, 3, 9, 8, 7};
    double test_doubles[] = {1.0, M_PI, M_PI*M_PI};
    char test_strings[] = "str1str2str3str4";
    tree.add_c_array(
        "test.int",
        test_ints,
        "int",
        6
    );
}
