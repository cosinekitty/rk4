#include <cstdio>
#include <cstring>

using test_func_t = int (*)();

struct Test
{
    const char *name;
    test_func_t func;
};


static int Happy();


static Test TestList[] =
{
    { "happy", Happy }
};


constexpr int NumTests = static_cast<int>(sizeof(TestList) / sizeof(TestList[0]));


int main(int argc, const char *argv[])
{
    if (argc == 2)
    {
        int rc;
        const char *verb = argv[1];
        if (!strcmp(verb, "all"))
        {
            for (int i = 0; i < NumTests; ++i)
            {
                printf("Running: %s\n", TestList[i].name);
                rc = TestList[i].func();
                if (rc)
                {
                    printf("FAIL: %s returned %d\n", TestList[i].name, rc);
                    return 1;
                }
            }
        }
        else
        {
            Test* test = nullptr;
            for (int i = 0; i < NumTests; ++i)
            {
                if (!strcmp(verb, TestList[i].name))
                {
                    test = &TestList[i];
                    break;
                }
            }
            if (test == nullptr)
            {
                printf("FAIL: Unknown test name '%s'\n", verb);
                return 1;
            }
            rc = test->func();
            if (rc)
            {
                printf("FAIL: %s returned %d\n", test->name, rc);
                return 1;
            }
        }
        printf("SUCCESS\n");
        return 0;
    }
    else
    {
        printf("unittest: Invalid command line arguments\n");
        return 1;
    }
}


static int Happy()
{
    printf("HAPPY!\n");
    return 0;
}
