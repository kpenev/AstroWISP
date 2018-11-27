#include "../SubPixPhot/CInterface.h"

#include "../IO/CommandLineConfig.h"

int main(int argc, char **argv)
{
    SubPixPhotConfiguration *config = create_subpixphot_configuration();

    update_subpixphot_configuration(
        config,
        ""
    );

    destroy_subpixphot_configuration(config);
}
