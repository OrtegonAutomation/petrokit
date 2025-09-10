from petrokit.nodal import nodal_analysis

def test_nodal_returns_valid_point():
    """
    El análisis nodal debe devolver un punto de operación válido.
    """
    q_op, pwf_op = nodal_analysis(
        p_res=3000,
        q_max=1200,
        well_depth=8000,
        rho=60,
        mu=1,
        d=2.992
    )

    assert q_op > 0
    assert pwf_op > 0
    assert pwf_op < 3000
    assert q_op < 1200 * 1.2


def test_nodal_consistency_with_qmax():
    """
    El caudal de operación no puede ser mayor que q_max,
    y si la VLP no es muy restrictiva, debería crecer con q_max.
    """
    q_op1, _ = nodal_analysis(
        p_res=3000,
        q_max=1000,
        well_depth=2000,  # pozo más somero → VLP menos restrictiva
        rho=60,
        mu=1,
        d=2.992
    )

    q_op2, _ = nodal_analysis(
        p_res=3000,
        q_max=2000,
        well_depth=2000,
        rho=60,
        mu=1,
        d=2.992
    )

    assert q_op1 > 0
    assert q_op2 > 0
    assert q_op1 <= 1000 * 1.2
    assert q_op2 <= 2000 * 1.2

    # en condiciones iguales de restricción, debería aumentar
    assert q_op2 > q_op1
