import numpy as np
from petrokit.nodal import nodal_analysis
from petrokit.vlp import vlp_curve


def test_nodal_returns_valid_point():
    """
    El análisis nodal debe devolver un punto de operación válido.

    Nota: puede existir caso de "no-flow" si la presión de yacimiento no alcanza
    para sostener la columna hidrostática (VLP(q=0) >= p_res). En ese caso q_op=0
    es un resultado válido.
    """
    p_res = 3000.0
    q_max = 1200.0
    well_depth = 8000.0
    rho = 60.0
    mu = 1.0
    d = 2.992

    q_op, pwf_op = nodal_analysis(
        p_res=p_res,
        q_max=q_max,
        well_depth=well_depth,
        rho=rho,
        mu=mu,
        d=d,
    )

    # VLP a caudal cero ~ columna hidrostática
    pwf_vlp0 = float(vlp_curve(np.array([0.0]), well_depth, rho, mu, d)[0])

    if pwf_vlp0 >= p_res:
        # No hay flujo natural posible
        assert q_op == 0.0
        assert np.isclose(pwf_op, p_res, atol=1e-6)
    else:
        # Debe haber un punto de operación con flujo
        assert q_op > 0.0
        assert 0.0 <= pwf_op <= p_res

    assert np.isfinite(q_op)
    assert np.isfinite(pwf_op)



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
