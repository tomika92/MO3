#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

void print(int n, double** tab)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << tab[i][j];
            cout << "\t\t";
        }
        cout << endl;
    }
}

void jacoby(const int n, double** tab, double* r, double eps)
{
    double rev[4];
    double wek2[4];
    double por[4];
    double wynik[4];
    double** tab2;
    double max;
    
    for (int i = 0; i < n; i++)
    {
        wynik[i] = r[i] / tab[i][i];
    }

    do
    {
        static int licznik = 0;
        for (int i = 0; i < n; i++)
        {
            rev[i] = wynik[i];
        }

        for (int j = 0; j < n; j++)
        {
            double sum = 0;
            for (int k = 0; k < n; k++)
            {
                if (j != k)
                {
                    sum = sum + tab[j][k] * rev[k];
                }
            }
            wynik[j] = (1.0 / tab[j][j]) * (r[j] - sum);
        }
        cout << "\nIteracja nr " << licznik + 1 << endl;
        cout << setprecision(6) << scientific;
        cout << "Wynik" << endl;
        for (int i = 0; i < n; i++)
        {
            cout << wynik[i] << "\t";
        }
        cout << endl;

        for (int i = 0; i < n; i++)
        {
            por[i] = abs(rev[i] - wynik[i]);
        }

        max = 0;
        for (int i = 0; i < n; i++)
        {
            if (por[i] > max)
            {
                max = por[i];
            }
        }
        cout << "Blad = " << max << endl;

        licznik++;
    } while (max>=eps);

    
    for (int i = 0; i < n; i++) //A*x
    {
        double s = 0;
        for (int j = 0; j < n; j++)
        {
            s += tab[i][j] * wynik[j];
        }
        por[i] = s;
    }

    for (int i = 0; i < n; i++)//|a*x-r|
    {
        wek2[i] = abs(por[i] - r[i]);
    }
   
    double max2 = 0;
    for (int i = 0; i < n; i++)
    {
        if (wek2[i] > max2)
        {
            max2 = wek2[i];
        }
    }

    cout << setprecision(6) << scientific;
    cout << "\n\nDokladnosc: " << max2 << endl;
}

void gauss_seidel(const int n, double** tab, double* r, double eps)
{
    double rev[4];
    double wynik[4];
    double por[4];
    double max;
    
    for (int i = 0; i < n; i++)
    {
        wynik[i] = r[i] / tab[i][i];
    }

    do
    {
        static int licznik = 0;
        for (int i = 0; i < n; i++)
        {
            rev[i] = wynik[i];
        }

        for (int i = 0; i < n; i++)
        {
            double sum1 = 0;
            double sum2 = 0;
            for (int j = 0; j < i; j++)
            {
                sum1 = sum1 + tab[i][j] * wynik[j];
            }
            for (int j = i + 1; j < n; j++)
            {
                sum2 = sum2 + tab[i][j] * rev[j];
            }
            wynik[i] = (1.0 / tab[i][i]) * (r[i] - sum1 - sum2);
        }
        cout << "\nIteracja nr " << licznik + 1 << endl;
        cout << "Wynik" << endl;
        for (int i = 0; i < n; i++)
        {
            cout << wynik[i] << "\t";
        }
        cout << endl;

        for (int i = 0; i < n; i++)
        {
            por[i] = abs(rev[i] - wynik[i]);
        }

        max = 0;
        for (int i = 0; i < n; i++)
        {
            if (por[i] > max)
            {
                max = por[i];
            }
        }
        cout << "Blad = " << max << endl;

        licznik++;
    } while (max >= eps);

    for (int i = 0; i < n; i++) //A*x
    {
        double s = 0;
        for (int j = 0; j < n; j++)
        {
            s += tab[i][j] * wynik[j];
        }
        por[i] = s;
    }

    for (int i = 0; i < n; i++)//|a*x-r|
    {
        rev[i] = abs(por[i] - r[i]);
    }

    double max2 = 0;
    for (int i = 0; i < n; i++)
    {
        if (rev[i] > max2)
        {
            max2 = rev[i];
        }
    }

    cout << setprecision(6) << scientific;
    cout << "\n\nDokladnosc: " << max2 << endl;
}

int main()
{
    const int n = 2;
    double eps = 0.0000001;
    double** tab;
    double r[n] = { 0.6804, 0.3872 };
    tab = new double* [n];
    ifstream file;
    file.open("dane.txt", ios::in);
    for (int i = 0; i < n; i++)
    {
        tab[i] = new double[n + 1];
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            file >> tab[i][j];
        }
    }
    file.close();

    cout << "\nMacierz pierwotna A" << endl;
    cout << setprecision(2) << fixed;
    print(n, tab);

    cout << "\nWektor r" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << r[i];
        cout << "\t";
    }
    cout << "\n\nMetoda Jacobiego" << endl;
    jacoby(n, tab, r, eps);
    
    cout << "\n\nMetoda Gaussa-Seidela" << endl;
    gauss_seidel(n, tab, r, eps);

    for (int i = 0; i < n; i++)
    {
        delete[] tab[i];
    }
}
