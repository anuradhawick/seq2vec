#pragma once
#include <iostream>
#include <iomanip>

using namespace std;

class ProgressDisplay
{
private:
    int total = 0;
    int progress = 0;
    int interval = 0;
    float percentage = 0.0;
public:
    ProgressDisplay(int total=0): total(total), interval(max(1, total/1000)){}

    void operator++(int)
    {
        operator++();
    }

    void operator++()
    {
        progress++;
        percentage = 100.0 * static_cast<float>(progress)/static_cast<float>(total);

        if (total % interval == 0)
        {
            print();
        }
    }

    void end()
    {
        progress = total;
        if (total > 0) {
            cout << "Completed " << fixed << setprecision(2) << 100.00 << "%       " << endl << flush;
        } else {
            cout << "Completed " << fixed << setprecision(2) << progress << "       " << endl << flush;
        }
    }

    void print()
    {
        if (total > 0) {
            cout << "Completed " << fixed << setprecision(2) << percentage << "%             \r" << flush;
        } else {
            cout << "Completed " << fixed << setprecision(2) << progress << "            \r" << flush;
        }
    }
};
