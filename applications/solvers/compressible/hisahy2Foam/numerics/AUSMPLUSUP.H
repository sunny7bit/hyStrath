
void AUSMPLUSUPFlux
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
    // const scalar alpha = 3.0/16.0;
    const scalar beta  = 1.0/8.0;
	const scalar Kp    = 0.25;
	const scalar Ku    = 0.75;
	const scalar sigma = 1.0;
    const scalar MaInf = 0.05;

    aByU = (U_pos + U_neg)/2.0;
    amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + 0.5*(cSf_pos + cSf_neg);
    // 输出中间速度计算
    // Pout << "Intermediate velocities:" << endl
    //      << "  aByU = " << aByU << endl
    //      << "  amaxSf (pre-update) = " << amaxSf << endl;
    //amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + c_f * magSf;
    //amaxSf = 0.5*mag((U_pos + U_neg) & Sf) + 0.5*(cSf_pos + cSf_neg);
    //amaxSf = max(mag(phiv_pos), mag(phiv_neg));
    // bounding variables
    //const scalar rhoMin = SMALL;

    // speed of sound, for left and right side, assuming perfect gas
    //const scalar c_pos  = Foam::sqrt(max(SMALL,kappaLeft *p_pos /max(rho_pos, rhoMin)));
    //const scalar c_neg = Foam::sqrt(max(SMALL,kappaRight*p_neg/max(rho_neg,rhoMin)));

    // DensityTotalEnergy
    const scalar rhoH_pos = rho_pos*(e_pos + 0.5 * magSqr(U_pos)) + p_pos;
    const scalar rhoH_neg = rho_neg*(e_neg + 0.5 * magSqr(U_neg)) + p_neg;

    // Step 1: decode left and right:

    // Compute conservative variables assuming perfect gas law

    // DensityVelocity
    //const vector rhoU_pos  = rho_pos *U_pos;
    //const vector rhoU_neg = rho_neg*U_neg;

    // compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    // normal vector
    const vector unit_normal = Sf / magSf;
    const scalar qLeft  = U_pos & unit_normal;
    const scalar qRight = U_neg & unit_normal;
    // 添加速度分量输出
    // Pout << "Velocity components:" << endl
    //      << "  qLeft = " << qLeft << ", qRight = " << qRight << endl;
    //方法一，界面统一声速用简单平均关系式，fa由计算获得
    // const scalar aTilde = 0.5*(c_pos+c_neg);
	// const scalar rhoTilde = 0.5*(rho_pos+rho_neg);
	// const scalar sqrMaDash = (sqr(qLeft)+sqr(qRight))/(2.0*sqr(aTilde));
	// const scalar sqrMaZero = min(1.0,max(sqrMaDash,sqr(MaInf)));
	// const scalar MaZero    = Foam::sqrt(sqrMaZero);
	// const scalar fa = MaZero*(2.0-MaZero);
	//方法二，总焓守恒下的界面统一声速，fa取定值1.0
    const scalar H_pos = rhoH_pos/rho_pos;
    const scalar H_neg = rhoH_neg/rho_neg;
    const scalar sqraStar_pos = 2.0*(gamma_pos-1)/(gamma_pos+1)*H_pos;
    const scalar sqraStar_neg = 2.0*(gamma_neg-1)/(gamma_neg+1)*H_neg;
    const scalar aHatLeft  = sqraStar_pos / max(Foam::sqrt(sqraStar_pos), qLeft);
    const scalar aHatRight = sqraStar_neg / max(Foam::sqrt(sqraStar_neg),-qRight);
    const scalar aTilde = min(aHatLeft, aHatRight);
	const scalar rhoTilde = 0.5*(rho_pos+rho_neg);
	const scalar sqrMaDash = (sqr(qLeft)+sqr(qRight))/(2.0*sqr(aTilde));
	// const scalar fa = 1.0;
    const scalar sqrMaZero = min(1.0,max(sqrMaDash,sqr(MaInf)));
	const scalar MaZero    = Foam::sqrt(sqrMaZero);
	const scalar fa = MaZero*(2.0-MaZero);
    // const scalar H_pos = rhoH_pos/rho_pos;
    // const scalar H_neg = rhoH_neg/rho_neg;
    // const scalar gamma = 0.5*(gamma_pos+gamma_neg);
    // const scalar sqraStar = 2.0*(gamma-1)/(gamma+1)*0.5*(H_pos+H_neg);
    // const scalar aStar = Foam::sqrt(sqraStar);
    // const scalar aHatLeft  = sqraStar / max(aStar, qLeft);
    // const scalar aHatRight = sqraStar / max(aStar,-qRight);
    // const scalar aTilde = min(aHatLeft, aHatRight);
	// const scalar rhoTilde = 0.5*(rho_pos+rho_neg);
	// const scalar sqrMaDash = (sqr(qLeft)+sqr(qRight))/(2.0*sqr(aTilde));
	// const scalar fa = 1.0;
    // Pout << "aTilde, qLeft, qRight, aHatLeft, aHatRight:" << endl
    //      << aTilde << ", " << qLeft << ", " << qRight << ", " << aHatLeft << ", " << aHatRight << endl;

	const scalar alpha = 3.0/16.0*(-4.0+5.0*sqr(fa));

    const scalar MaRelLeft  = qLeft / aTilde;
    const scalar MaRelRight = qRight / aTilde;

    const scalar magMaRelLeft  = mag(MaRelLeft);
    const scalar magMaRelRight = mag(MaRelRight);

    const scalar Ma1PlusLeft   = 0.5*(MaRelLeft +magMaRelLeft );
    const scalar Ma1MinusRight = 0.5*(MaRelRight-magMaRelRight);

    const scalar Ma2PlusLeft   =  0.25*sqr(MaRelLeft +1.0);
    const scalar Ma2PlusRight  =  0.25*sqr(MaRelRight+1.0);
    const scalar Ma2MinusLeft  = -0.25*sqr(MaRelLeft -1.0);
    const scalar Ma2MinusRight = -0.25*sqr(MaRelRight-1.0);

    const scalar Ma4BetaPlusLeft   = ((magMaRelLeft  >= 1.0) ? Ma1PlusLeft   : (Ma2PlusLeft  *(1.0-16.0*beta*Ma2MinusLeft)));
    const scalar Ma4BetaMinusRight = ((magMaRelRight >= 1.0) ? Ma1MinusRight : (Ma2MinusRight*(1.0+16.0*beta*Ma2PlusRight)));
	
	const scalar Mp = -Kp/fa*max(1.0-sigma*sqrMaDash,0.0)*(p_neg-p_pos)/(rhoTilde*sqr(aTilde));

    const scalar P5alphaPlusLeft   = ((magMaRelLeft  >= 1.0) ?
        (Ma1PlusLeft/MaRelLeft)    : (Ma2PlusLeft  *(( 2.0-MaRelLeft) -16.0*alpha*MaRelLeft *Ma2MinusLeft )));
    const scalar P5alphaMinusRight = ((magMaRelRight >= 1.0) ?
        (Ma1MinusRight/MaRelRight) : (Ma2MinusRight*((-2.0-MaRelRight)+16.0*alpha*MaRelRight*Ma2PlusRight)));
	
	const scalar pU = -Ku*P5alphaPlusLeft*P5alphaMinusRight*(rho_pos+rho_neg)*(fa*aTilde)*(qRight-qLeft);
	    
    const scalar MaRelTilde = Ma4BetaPlusLeft + Ma4BetaMinusRight + Mp;
    const scalar pTilde = p_pos*P5alphaPlusLeft + p_neg*P5alphaMinusRight + pU;
    // 输出最终组合值
    // Pout << "Combined values:" << endl
    //      << "  MaRelTilde = " << MaRelTilde << ", pTilde = " << pTilde << endl;

    const scalar URelTilde = MaRelTilde*aTilde;
    const scalar magURelTilde = mag(MaRelTilde)*aTilde;
    // 关键检查点4：速度计算
    // Pout << "Interface velocity:" << endl
    //      << "  URelTilde = " << URelTilde << ", magURelTilde = " << magURelTilde << endl;

    // There is a typo in Luo et. al, J. Comp. Physics 194 (2004), Chap 4.2 Eq. 4.8
    // refer to the origial Paper from Liou, J. Comp. Physics 129 (1996), Chap4, Eq. 42
    rhoFlux  = 0.5*(URelTilde*(rho_pos +rho_neg) -magURelTilde*(rho_neg -rho_pos))*magSf;
    rhoUFlux = 0.5*(URelTilde*(rhoU_pos+rhoU_neg)-magURelTilde*(rhoU_neg-rhoU_pos))*magSf+pTilde*Sf;
    rhoEFlux = 0.5*(URelTilde*(rhoH_pos+rhoH_neg)-magURelTilde*(rhoH_neg-rhoH_pos))*magSf;
    // 最终通量输出
    // Pout << "Final fluxes:" << endl
    //      << "  rhoFlux = " << rhoFlux << endl
    //      << "  rhoUFlux = " << rhoUFlux << endl
    //      << "  rhoEFlux = " << rhoEFlux << endl
    //      << "===========================================" << endl << endl;
    
    // ---
    if(downgradeSingleT)
    {
        rhoEvkFlux.clear();
    }
    else if(downgradeSingleTv)
    {
        rhoEvkFlux[0] = 0.5*(URelTilde*(rho_pos*evk_pos[0] +rho_neg*evk_neg[0]) -magURelTilde*(rho_neg*evk_neg[0] -rho_pos*evk_pos[0]))*magSf;
    }
    else
    {
        forAll(rhoEvkFlux, speciei)
        { 
            rhoEvkFlux[speciei] = 0.5*(URelTilde*(pD_pos[speciei]*evk_pos[speciei] +pD_neg[speciei]*evk_neg[speciei]) -magURelTilde*(pD_neg[speciei]*evk_neg[speciei] -pD_pos[speciei]*evk_pos[speciei]))*magSf;
        }
    }
}

