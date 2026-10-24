#include <cstdio>
#include <cstring>
#include "SGP4.h"

struct Case { const char* id; bool shouldAccept; int expectSatnum; };

int main()
{
    Case cases[] = {
        { "00005", true,   5 },
        { "25544", true,   25544 },
        { "E8493", true,   148493 },
        { "O1234", false, -1 },
        { "I0000", false, -1 },
        { "AA000", false, -1 },
        { "#1234", false, -1 },
    };

    const char line1t[] = "1 00005U 58002D   26189.70990935 -.00000023  00000+0 -70517-5 0  9996";
    const char line2t[] = "2 00005  34.2417 341.8745 1487004  19.9191 345.3718 11.62373363189308";

    int failures = 0;
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++)
    {
        char line1[130], line2[130];
        strcpy(line1, line1t);
        strcpy(line2, line2t);
        memcpy(line1 + 2, cases[i].id, 5);
        memcpy(line2 + 2, cases[i].id, 5);

        double startmfe, stopmfe, deltamin;
        elsetrec satrec;
        memset(&satrec, 0, sizeof(satrec));
        SGP4Funcs::twoline2rv(line1, line2, 'c', 'm', 'i', wgs72, startmfe, stopmfe, deltamin, satrec);

        bool passed = cases[i].shouldAccept
            ? (satrec.error == 0 && satrec.satnum == cases[i].expectSatnum)
            : (satrec.error == 7 && satrec.satnum == -1);

        printf("%-6s expect %-7s -> satnum=%-7d error=%-2d  %s\n",
            cases[i].id, cases[i].shouldAccept ? "accept" : "refuse",
            satrec.satnum, satrec.error, passed ? "PASS" : "FAIL");

        if (!passed)
            failures++;
    }

    printf("\n%d case(s) failed\n", failures);
    return failures != 0;
}
