import numpy as np
import pytest

from petrokit.nodal import nodal_analysis


def test_nodal_default_still_works():
    q_op, pwf_op = nodal_analysis(
        p_res=3000.0,
        q_max=1200.0,
        well_depth=8000.0,
        rho=60.0,
        mu=1.0,
        d=2.992,
        npts=50,
    )
    assert np.isfinite(q_op)
    assert np.isfinite(pwf_op)


def test_nodal_invalid_model_raises():
    with pytest.raises(ValueError):
        nodal_analysis(
            p_res=3000.0,
            q_max=1200.0,
            well_depth=8000.0,
            rho=60.0,
            mu=1.0,
            d=2.992,
            npts=50,
            ipr_model="no_existe",
        )
