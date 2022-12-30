
void GenerateJac()
{
    int i, j, k, l, m;
    int nElm, nonzeros_B;
    int Jac_SP, Jac;

    /* Each nonzero entry of B now counts its rank */
    nonzeros_B = 0;
    for (i = 0; i < EqnNr; i++) {
        for (j = 0; j < SpcNr; j++) {
            if (structB[i][j] != 0)
            {
                nonzeros_B++;
                structB[i][j] = nonzeros_B;
            }
        }
    }

    if ((useLang == C_LANG) || (useLang == F77_LANG) || (useLang == F90_LANG))
    {
        NewLines(1);
        WriteComment("Local variables");
        /* DeclareConstant( NTMPB,   ascii( nonzeros_B ) ); */
        varTable[NTMPB]->value = nonzeros_B;
        Declare(BV);
    }

    for (i = 0; i < EqnNr; i++)
    {
        for (j = 0; j < VarNr; j++)
        {
            if (Stoich_Left[j][i] != 0)
            {
                prod = Mul(RConst(i), Const(Stoich_Left[j][i]));
                for (l = 0; l < VarNr; l++)
                {
                    m = (int)Stoich_Left[l][i] - (l == j);
                    for (k = 1; k <= m; k++)
                        prod = Mul(prod, Elm(V, l));
                }
                for (; l < SpcNr; l++)
                    for (k = 1; k <= (int)Stoich_Left[l][i]; k++)
                        prod = Mul(prod, Elm(F, l - VarNr));
                /* Comment the B */
                WriteComment("B(%d) = dA(%d)/dV(%d)", Index(structB[i][j] - 1), Index(i), Index(j));
                Assign(Elm(BV, structB[i][j] - 1), prod);
            }
        }
    }

    nElm = 0;
    NewLines(1);
    WriteComment("Construct the Jacobian terms from B's");

    for (i = 0; i < VarNr; i++)
    {
        for (j = 0; j < VarNr; j++)
        {
            if (structJ[i][j])
            {
                sum = Const(0);
                for (k = 0; k < EqnNr; k++)
                {
                    if (Stoich[i][k] * structB[k][j] != 0)
                        sum = Add(sum, Mul(Const(Stoich[i][k]), Elm(BV, structB[k][j] - 1)));
                }
                Assign(Elm(JV, i, j), sum);
            }
        }
    }

}
