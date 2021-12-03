%:- consult(KB).

%goal(S).

isMove(X):-
    X = up;
    X = down;
    X = left;
    X = right;
    X = carry;
    X = drop.

result(M, s0):-
    isMove(M).

result(M0, result(M1, X)):-
    isMove(M0),
    isMove(M1),
    result(M1, X).