"""Checks for competition, overlap, time recurrence and seam handling."""
import numpy as np
from .core import graph,intervals,temporal_recurrence,shifted_targets,unique_count


def event(peak,polarity=1):
    return dict(start=peak-1,peak=peak,centroid=peak,end=peak+1,polarity=polarity)


def test_competing_source_cannot_be_unique():
    s=[event(10),event(15)];g=[event(20)]
    edges,sd,td=graph(s,g,60)
    assert edges==[(0,0),(1,0)] and list(sd)==[1,1] and list(td)==[2]
    assert unique_count(s,g,60)==0


def test_all_competing_targets_retained():
    edges,sd,td=graph([event(10)],[event(15),event(20)],60)
    assert len(edges)==2 and sd[0]==2
    assert unique_count([event(10)],[event(15),event(20)],60)==0


def test_direction_polarity_and_endpoints():
    s=[event(10)];g=[event(9),event(10),event(70),event(71),event(20,-1)]
    assert graph(s,g,60)[0]==[(0,1),(0,2)]
    assert graph(s,g,60,False)[0]==[(0,1),(0,2),(0,4)]


def test_transitive_overlap_not_independent_samples():
    assert intervals([event(10),event(60),event(110),event(200)],60)==[0,0,0,1]


def test_recurrence_never_double_assigns():
    ref=[event(10),event(12),event(30,-1)]
    assert len(temporal_recurrence(ref,[event(11)],5))==1
    assert temporal_recurrence([event(10,-1)],[event(10)],5)==set()


def test_shift_preserves_lag_shape_and_drops_seam_lobes():
    result=shifted_targets([event(10),event(89),event(90)],10,100)
    assert len(result)==1 and result[0]['peak']==20
    assert result[0]['end']-result[0]['start']==2


def test_replay_calibration_interval_matches_frozen_kernel():
    import json
    from pathlib import Path
    from ..response_broadening.core import KernelDesign,changes,detect_events
    from .run import events_at,RESEARCH,PRIOR
    cfg=json.loads((RESEARCH/'config.json').read_text())
    with np.load(PRIOR/'baseline/data/6B03205_q2_processed.npz') as z:
        expected=detect_events(changes(z['oxygen_flow'],cfg),cfg,KernelDesign(changes(z['oxygen_flow'],cfg),cfg).train)
        actual,_=events_at(z,'co',cfg,31,9)
        assert len(actual)==len(expected)==111
        assert [r['peak'] for r in actual]==[r['peak'] for r in expected]


def test_recurrence_maximizes_valid_pairs_before_distance():
    # Raw distance minimization can discard a valid pair in favor of two
    # out-of-tolerance assignments. Invalid edges must carry a fixed penalty.
    assert temporal_recurrence([event(10),event(20)],[event(14),event(4)],5)=={0}
