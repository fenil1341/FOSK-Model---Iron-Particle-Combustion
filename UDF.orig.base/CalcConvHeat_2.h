#include "udf.h"
#include "dpm.h"

real calc_conv_heat(Tracked_Particle *p)
{
    Thread *t = TP_CELL_THREAD(p);
    cphase_state_t *c = &(p->cphase[0]);
    //static int counter = 0;

    real p_ID = TP_ID(p);
    real Tp = TP_T(p);            /* particle temperature */
    real mp = TP_MASS(p);         /* particle mass */
    real Dp = DPM_DIAM_FROM_VOL(mp / TP_RHO(p));    /* particle diameter */
    real Ap = DPM_AREA(Dp);                        /* particle surface */
    real Pr = c->sHeat * c->mu / c->tCond;         /* Prandtl number */
    real Re = p->Re;
    real Nu = 2.0 + 0.6 * sqrt(p->Re) * pow(Pr, 1./3.);  /* Nusselt number */
    real h = Nu * c->tCond / Dp;                 /* Heat transfer coefficient*/
    real Tg = c->temp;
    real dh_dt = h * (c->temp - Tp) * Ap;        /* heat source term*/
    
    real dt = TP_DT(p);
    real cp_p = DPM_SPECIFIC_HEAT(p, Tp);
    real dT_pred = dh_dt * dt / (mp * cp_p);
    real limiter = 0.0;
    real dT_pred_2 = 0.0;
    

    if (fabs(dT_pred) > fabs(Tg - Tp))
    {
        dT_pred = Tg - Tp;

        dh_dt = mp * cp_p * dT_pred / dt;

        limiter = 1.0;

        dT_pred_2 = dh_dt * dt / (mp * cp_p);

    }


    if (fabs(Tg - Tp) < 1e-1)
    {
        dh_dt = 0.0;

        limiter += 2.0;
    }
    /*
    
    counter++;
    if (counter % 50000 == 0)
    {
        Message("dT_pred=%.5e | dT_pred_2=%.5e | Pr=%.5e  | h=%.5e   | Dp=%.5e | mp=%.5e | Ap=%.5e |Re=%.5e | Nu=%.5e | Tg=%.5e | Tp=%.5e | cp_p=%.5e | dh_dt=%.5e | dt=%.5e | limiter=%.1f \n",
            dT_pred, dT_pred_2, Pr, h, Dp, mp, Ap, Re, Nu, Tg, Tp, cp_p, dh_dt, dt, limiter);
    }
    */

    
    return dh_dt;
}