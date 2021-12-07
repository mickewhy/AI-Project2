% Importing the knowledgebase.
:- include('KB.pl').

% __________MAIN__________

% Returns true if given a goal state or returns the goal state if given a variable.
%goal(result(M,X)):-
%    result(M,X),


% Base successor state, is true if given a valid move M and the initial state s0.
result(M, s0):-
    isMove(M).

% Recursive successor state, is true iff M0, M1 are valid moves AND result(M1, X) is true.
result(M0, result(M1, X)):-
    isMove(M0),
    isMove(M1),
    result(M1, X).
    

% __________HELPERS__________

% Checks if X is a valid move.
% ?- isMove(up). true ; isMove(north). false
isMove(X):-
    X = up;
    X = down;
    X = left;
    X = right;
    X = carry;
    X = drop.

% Returns move that gets NX,NY closer to HX,HY. Returns false if coordinates are equal.
% ?- moveTo(0,0,1,0,A). A = down
moveTo(NX,NY,HX,HY,A):-
    (HX>NX , A = down);
    (HX<NX , A = up);
    (HY>NY, A = right);
    (HY<NY, A = left).

% If NeoXY = HostageXY, state = s0.
getMovesTo(X,Y,X,Y,s0).

% Otherwise, state = result(action, newState).
% ?- getMovesTo(0,0,1,1,M). M = result(down, result(right, s0)) 
getMovesTo(NX,NY,HX,HY,M):-
    moveTo(NX,NY,HX,HY,A),
    ((A = up, NEW is NX-1, getMovesTo(NEW,NY,HX,HY,M1)) ; (A = down, NEW is NX+1, getMovesTo(NEW,NY,HX,HY,M1)) ; (A = right, NEW is NY+1, getMovesTo(NX,NEW,HX,HY,M1)) ; (A = left, NEW is NY-1, getMovesTo(NX,NEW,HX,HY,M1))),
    M = result(A,M1).

% Base case, if given state is s0 then the hostages in s0 = the length of hostages_loc.
% ?- calcHostages(s0,X). X = 2
calcHostages(s0, X):-
    hostages_loc(L),
    length(L,X).

% ?- calcHostages(result(drop,s0), X). X = 1
calcHostages(result(M,X), H):-
    result(M,X),
    calcHostages(X, H1),
    ((M = drop, H is H1-1);(H is H1)).

% Base case, if given state is s0 then the capacity X in s0 = capacity(X).
% ?- calcCapacity(s0,X). X = 1
calcCapacity(s0, C):-
    capacity(C).

% calcCapacity(result(carry,s0),X). X = 0
calcCapacity(result(M,X), C):-
    result(M,X),
    calcCapacity(X, C1),
    ((M = drop, C is C1+1);(M = carry, C is C1-1);(C is C1)).

% Base case, if given state is s0 then Neo's position is neo_loc
% ?- calcNeo(s0,X,Y). X = Y, Y = 0
calcNeo(s0,NX,NY):-
    neo_loc(NX,NY).

% ?- calcNeo(result(down,s0),X,Y). X = 1, Y = 0
% ?- calcNeo(result(right,s0),X,Y). X = 0, Y = 1
% ?- calcNeo(result(right,result(down,s0)),X,Y). X = Y, Y = 1
calcNeo(result(M,X), NX, NY):-
    result(M,X),
    calcNeo(X,NX1,NY1),
    ((M = up, NX is NX1-1, NY is NY1);(M = down, NX is NX1+1, NY is NY1);(M = right, NX is NX1, NY is NY1+1);(M = left, NX is NX1, NY is NY1-1)).

% Gets hostage depending on the number of hostages H.
% If H = 2, like in the KB, this will bind HXY to the hostage at location(H1 = H-1 = 1), i.e. the second hostage.
% If H = 1, i.e. we only have 1 hostage, this will bind HXY to the hostage at location(H1 = H-1 = 0), i.e. the first hostage.
% ?- getHostage(s0,H). H = [1, 2]
% ?- getHostage(result(drop,s0),H). H = [1, 1]
getHostage(S,HXY):-
    calcHostages(S,H),
    hostages_loc(L),
    H1 is H-1,
    nth0(H1, L, HXY).

% Base case, checks if all of the hostages are saved.
% X can be s0 with empty list of hostages or recursive result() state that leads to all hostages being saved.
matrix(X):-
    calcHostages(X,0).

matrix(R, CARRY, OUT):-
    R = result(M,X),
    calcHostages(R, H),
    calcCapacity(R, C),
    calcNeo(R, NX, NY),
    booth(BX, BY),
    getHostage(R,[HX, HY]),
    getMovesTo(NX, NY, HX, HY, MOVESH),
    CARRY = result(carry,MOVESH),
    calcNeo(MOVESH, NX1, NY1),
    getMovesTo(NX1, NY1, BX, BY, MOVESB),
    OUT = result(drop,MOVESB).

