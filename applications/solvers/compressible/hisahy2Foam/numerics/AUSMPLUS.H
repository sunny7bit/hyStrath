
void AUSMPLUSFlux
(
    const bool& downgradeSingleT,
    const bool& downgradeSingleTv,
    const vector& Sf,
    const scalar& magSf,
    const scalar& c_pos,
    const scalar& c_neg,
    const scalar& gamma_pos,
    const scalar& gamma_neg,
    const scalar& cSf_pos,
    const scalar& cSf_neg,
    const scalar& rho_pos,
    const scalar& rho_neg,
    const PtrList<scalar>& pD_pos,
    const PtrList<scalar>& pD_neg,
    const vector& U_pos,
    const vector& U_neg,
    const vector& rhoU_pos,
    const vector& rhoU_neg,
    const scalar& p_pos,
    const scalar& p_neg,
    const PtrList<scalar>& evk_pos,
    const PtrList<scalar>& evk_neg,
    const scalar& e_pos,
    const scalar& e_neg,
    scalar& amaxSf,
    vector& aByU,
    scalar& rhoFlux,
    vector& rhoUFlux,
    PtrList<scalar>& rhoEvkFlux,
    scalar& rhoEFlux
)
{
    // 添加调试输出头
    // Pout << "========== AUSMPLUSFlux Debug Pout ==========" << endl;

    // 输出输入参数关键值
    // Pout << "Inputs:" << endl
    //      << "  c_pos = " << c_pos << ", c_neg = " << c_neg << endl
    //      << "  rho_pos = " << rho_pos << ", rho_neg = " << rho_neg << endl
    //      << "  p_pos = " << p_pos << ", p_neg = " << p_neg << endl
    //      << "  U_pos = " << U_pos << ", U_neg = " << U_neg << endl;

    const scalar alpha = 3.0/16.0;
    const scalar beta = 1.0/8.0;

    aByU = (U_pos + U_neg)/2.0;
    amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + 0.5*(cSf_pos + cSf_neg);
    // 输出中间速度计算
    // Pout << "Intermediate velocities:" << endl
    //      << "  aByU = " << aByU << endl
    //      << "  amaxSf (pre-update) = " << amaxSf << endl;

    const scalar rhoH_pos = rho_pos*(e_pos + 0.5 * magSqr(U_pos)) + p_pos;
    const scalar rhoH_neg = rho_neg*(e_neg + 0.5 * magSqr(U_neg)) + p_neg;

    const vector unit_normal = Sf / magSf;
    const scalar u_l = U_pos & unit_normal;
    const scalar u_r = U_neg & unit_normal;

    // 添加速度分量输出
    // Pout << "Velocity components:" << endl
    //      << "  u_l = " << u_l << ", u_r = " << u_r << endl;

    const scalar c_f = 0.5*(c_pos + c_neg);
    // const scalar c_f = Foam::sqrt(c_pos*c_neg);

    // 关键检查点1：声速计算
    // if (c_pos <= 0 || c_neg <= 0) {
    //     Pout << "ERROR: Invalid speed of sound! c_pos = " << c_pos 
    //          << ", c_neg = " << c_neg << endl;
    // }
    // Pout << "Computed c_f = " << c_f << endl;

    //amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + c_f * magSf;
    //amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + 0.5*(cSf_pos + cSf_neg);
    //amaxSf = max(mag(phiv_pos), mag(phiv_neg));

    // const scalar Ma_l = u_l/c_f;
    // const scalar Ma_r = u_r/c_f;
    const scalar Ma_l = u_l/c_f;
    const scalar Ma_r = u_r/c_f;

    // 关键检查点2：马赫数
    // Pout << "Mach numbers:" << endl
    //      << "  Ma_l = " << Ma_l << ", Ma_r = " << Ma_r << endl;
    // if (Foam::mag(Ma_l) > 10 || Foam::mag(Ma_r) > 10) {
    //     Pout << "WARNING: Extreme Mach number detected!" << endl;
    // }

    const scalar Ma_lmag = mag(Ma_l);
    const scalar Ma_rmag = mag(Ma_r);

    const scalar Ma_pl = ((Ma_lmag >= 1) ? (0.5*(Ma_l + Ma_lmag)) : (0.25*sqr(Ma_l + 1) + beta*sqr(sqr(Ma_l) - 1)));
    const scalar Ma_nr = ((Ma_rmag >= 1) ? (0.5*(Ma_r - Ma_rmag)) : (-0.25*sqr(Ma_r - 1) - beta*sqr(sqr(Ma_r) - 1)));

    // 输出分裂马赫数
    // Pout << "Split Mach numbers:" << endl
    //      << "  Ma_pl = " << Ma_pl << ", Ma_nr = " << Ma_nr << endl;

    const scalar p_pl = ((Ma_lmag >= 1) ? (0.5*(Ma_l + Ma_lmag)/Ma_l) : (0.25*sqr(Ma_l + 1)*(2 - Ma_l) + alpha*Ma_l*sqr(sqr(Ma_l) - 1)));
    const scalar p_nr = ((Ma_rmag >= 1) ? (0.5*(Ma_r - Ma_rmag)/Ma_r) : (0.25*sqr(Ma_r - 1)*(2 + Ma_r) - alpha*Ma_r*sqr(sqr(Ma_r) - 1)));

    // 关键检查点3：压力分裂项
    // Pout << "Pressure split terms:" << endl
    //      << "  p_pl = " << p_pl << ", p_nr = " << p_nr << endl;
    // if (Foam::mag(p_pl) > 1e6 || Foam::mag(p_nr) > 1e6) {
    //     Pout << "WARNING: Abnormal pressure split detected!" << endl;
    // }

    const scalar Ma_f = Ma_pl + Ma_nr;
    const scalar p_f = p_pl*p_pos + p_nr*p_neg;

    // 输出最终组合值
    // Pout << "Combined values:" << endl
    //      << "  Ma_f = " << Ma_f << ", p_f = " << p_f << endl;

    const scalar u_f = c_f * Ma_f;
    const scalar u_fmag = c_f * mag(Ma_f);

    // 关键检查点4：速度计算
    // Pout << "Interface velocity:" << endl
    //      << "  u_f = " << u_f << ", u_fmag = " << u_fmag << endl;
    // if (Foam::mag(u_f) > 1e6 || u_fmag < 0) {
    //     Pout << "ERROR: Invalid interface velocity!" << endl;
    // }

    rhoFlux = 0.5*(u_f*(rho_pos + rho_neg) - u_fmag*(rho_neg - rho_pos))*magSf;
    rhoUFlux = 0.5*(u_f*(rhoU_pos + rhoU_neg) - u_fmag*(rhoU_neg - rhoU_pos))*magSf + p_f*Sf;
    rhoEFlux = 0.5*(u_f*(rhoH_pos + rhoH_neg) - u_fmag*(rhoH_neg - rhoH_pos))*magSf;
    // 最终通量输出
    // Pout << "Final fluxes:" << endl
    //      << "  rhoFlux = " << rhoFlux << endl
    //      << "  rhoUFlux = " << rhoUFlux << endl
    //      << "  rhoEFlux = " << rhoEFlux << endl
    //      << "===========================================" << endl << endl;
    // // ---
    if(downgradeSingleT)
    {
        rhoEvkFlux.clear();
    }
    else if(downgradeSingleTv)
    {
        rhoEvkFlux[0] = 0.5*(u_f*(rho_pos*evk_pos[0] +rho_neg*evk_neg[0]) -u_fmag*(rho_neg*evk_neg[0] -rho_pos*evk_pos[0]))*magSf;
    }
    else
    {
        forAll(rhoEvkFlux, speciei)
        { 
            rhoEvkFlux[speciei] = 0.5*(u_f*(pD_pos[speciei]*evk_pos[speciei] +pD_neg[speciei]*evk_neg[speciei]) -u_fmag*(pD_neg[speciei]*evk_neg[speciei] -pD_pos[speciei]*evk_pos[speciei]))*magSf;
        }
    }
}

