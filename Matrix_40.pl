% Importing the knowledgebase.
:- include('KB.pl').

% __________MAIN__________

% Returns true if given a goal state or returns the goal state if given a variable.
% Base case, checks if all of the hostages are saved, i.e. hostages = 0.
% X can be s0 with empty list of hostages or recursive result() state that leads to all hostages being saved.
% ?- goal(result(drop, result(up, result(carry,result(down,result(drop, result(up, result(right, result(carry,result(down, result(right, s0))))))))))). true
% ?- goal(result(up, result(carry, result(down,result(drop, result(up, result(right, result(carry,result(down, result(right, s0)))))))))). false
% ?- goal(S). S = result(drop, result(up, result(right, result(carry, result(down, result(left, result(drop, result(up, result(carry, result(down, result(right, result(right, s0))))))))))))
goal(X):-
    % there are maximum two hostages.
    hostages_loc(L), length(L,V), @=<(V,2),

    % the grid’s size is either 3 × 3 or 4 × 4.
    grid(M,M),((M = 3);(M = 4)),

    % calculate goal
    (\+var(X), calcHostages(X,0), assertLocations(X)) ; (var(X), goalRec(s0,X)).

% Base successor state, is true if given a valid move M and the initial state s0.
% ?- result(M,s0).
% M = up ;
% M = down ;
% M = left ;
% M = right ;
% M = carry ;
% M = drop.
result(M, s0):-
    isMove(M).

% Recursive successor state, is true iff M0, M1 are valid moves AND result(M1, X) is true.
% ?- result(up, result(down, s0)). true
% ?- result(up, result(south, s0)). false
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
% ?- moveTo(0,0,0,0,A). false
moveTo(NX,NY,HX,HY,A):-
    (HX>NX , A = down);
    (HX<NX , A = up);
    (HY>NY, A = right);
    (HY<NY, A = left).

% Gets the state of results that moves the first coordinate set to the second coordinate set
% If NeoXY = HostageXY, state = s0.
getMovesTo(X,Y,X,Y,s0).

% Otherwise, state = result(action, newState).
% ?- getMovesTo(0,0,1,1,M).
% M = result(down, result(right, s0)) ;
% M = result(right, result(down, s0)) ;
% false.
getMovesTo(NX,NY,HX,HY,M):-
    moveTo(NX,NY,HX,HY,A),
    ((A = up, NEW is NX-1, getMovesTo(NEW,NY,HX,HY,M1)) ; (A = down, NEW is NX+1, getMovesTo(NEW,NY,HX,HY,M1)) ; (A = right, NEW is NY+1, getMovesTo(NX,NEW,HX,HY,M1)) ; (A = left, NEW is NY-1, getMovesTo(NX,NEW,HX,HY,M1))),
    M = result(A,M1).

% Caculates the number of hostages in a given state
% Base case, if the given state is s0 then the hostages in s0 = the length of hostages_loc.
% ?- calcHostages(s0,X). X = 2
calcHostages(s0, X):-
    hostages_loc(L),
    length(L,X).

% Caculates the capacity in a given state
% Recursively calculates the number of hostages by subtracting 1 with each drop
% ?- calcHostages(result(drop,s0), X). X = 1
calcHostages(result(M,X), H):-
    result(M,X),
    calcHostages(X, H1),
    ((M = drop, H is H1-1);(H is H1)).

% Base case, if given state is s0 then the capacity X in s0 = capacity(X).
% ?- calcCapacity(s0,X). X = 1
calcCapacity(s0, C):-
    capacity(C).

% Recursively calculates the capacity by subtracting 1 with each carry and adding 1 with each drop
% calcCapacity(result(carry,s0),X). X = 0
calcCapacity(result(M,X), C):-
    result(M,X),
    calcCapacity(X, C1),
    ((M = drop, C is C1+1);(M = carry, C is C1-1);(C is C1)).

% Caculates Neo's coordinates in a given state
% Base case, if the given state is s0 then Neo's position is neo_loc
% ?- calcNeo(s0,X,Y). X = Y, Y = 0
calcNeo(s0,NX,NY):-
    neo_loc(NX,NY).

% Recursively calculates Neo's coordinates by calculating the lesser state's coordinates, then changing them according to the current move
% ?- calcNeo(result(down,s0),X,Y). X = 1, Y = 0
% ?- calcNeo(result(right,s0),X,Y). X = 0, Y = 1
% ?- calcNeo(result(right,result(down,s0)),X,Y). X = Y, Y = 1
calcNeo(result(M,X), NX, NY):-
    result(M,X),
    calcNeo(X,NX1,NY1),
    ((M = up, NX is NX1-1, NY is NY1);(M = down, NX is NX1+1, NY is NY1);(M = right, NX is NX1, NY is NY1+1);(M = left, NX is NX1, NY is NY1-1);((M=carry;M=drop), NX is NX1, NY is NY1)).

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

% Merges two states and returns them in R.
% ?- merge(s0, result(drop,s0), R). R = result(drop, s0).
merge(s0, M1, R):-
    R = M1.

% If the first state isn't s0, R = result(Mn,M1) where Mn = all of the moves originally in the first state.
% ?- merge(result(left,result(right,s0)),result(down,result(up,s0)),R). R = result(left, result(right, result(down, result(up, s0))))
merge(result(M,X),M1,R):-
    (X = s0, R = result(M,M1)) ; (merge(X, M1, R1), R = result(M,R1)).

% Drops one hostage from state R at the telephone booth and returns new state R1.
% ?- drop1(s0,R1). R1 = result(drop, result(up, result(carry, result(down, result(right, result(right, s0))))))
drop1(R,R1):-
    (R = result(_,_);R = s0),
    %calcHostages(R, H),
    %calcCapacity(R, C),
    calcNeo(R, NX, NY),
    booth(BX, BY),
    getHostage(R,[HX, HY]),
    getMovesTo(NX, NY, HX, HY, MOVESH),
    CARRY = result(carry,MOVESH),
    getMovesTo(HX, HY, BX, BY, MOVESB),
    DROP = result(drop,MOVESB),
    merge(DROP, CARRY, R1).

% Recursively calls drop1 as long as the hostage count isn't equal to zero.
% If the count reaches zero, return the merged states retrieved from drop1. Otherwise, call goalRec with those merged states and merge the result.
% ?- goalRec(s0, X). X = result(drop, result(up, result(right, result(carry, result(down, result(left, result(drop, result(up, result(carry, result(down, result(right, result(right, s0))))))))))))
goalRec(R,RET):-
    drop1(R,RTEMP),
    merge(RTEMP, R, R1),
    calcHostages(R1, H),
    ((H = 0, RET = R1) ; goalRec(R1,RET)).

% Checks if the drop and carry locations are correct in the given state.
% Base case, s0 has no drops or carries therefore assertLocations(s0) is a fact.
assertLocations(s0).

% Recursively calls lesser states until it reaches s0. With each state, check the move M.
% Do nothing if M is a move, i.e. M = up/down/left/right.
% If M is drop, calculate the current position for Neo and make sure he is standing at the telephone booth.
% If M is carry, calculate the current position for Neo and make sure he is standing at a hostage's location.
% ?- assertLocations(result(drop,s0)). false
% ?- assertLocations(result(drop, result(up, result(carry, result(down, result(right, result(right, s0))))))). true
% ?- assertLocations(result(carry,s0)). false
% ?- assertLocations(result(carry,result(down, result(right, s0)))). true
assertLocations(result(M,X)):-
    ( (M = drop, calcNeo(result(M,X), NX, NY), booth(NX,NY)) ; (M = carry, calcNeo(result(M,X), NX, NY), hostages_loc(L), member([NX,NY], L)) ; (M = up) ; (M = down) ; (M = left) ; (M = right) ),
    assertLocations(X).