graph [
  directed 1
  state ""
  meta ""
  mutype ""
  node [
    id 0
    label "0"
    state "[{a,a}, [b,b]]"
    mutype [
      hetA 2
      hetB 2
    ]
    meta "source"
      pos 207.0
      pos 306.0
  ]
  node [
    id 1
    label "1"
    state "[{}, [a,a,b,b]]"
    mutype [
      hetA 2
      hetB 2
    ]
    meta ""
      pos 135.0
      pos 234.0
  ]
  node [
    id 2
    label "2"
    state "[{aa}, [b,b]]"
    mutype [
      fixed 1
      hetB 2
    ]
    meta ""
      pos 207.0
      pos 234.0
  ]
  node [
    id 3
    label "3"
    state "[{a,a}, [bb]]"
    mutype [
      hetA 2
      fixed 1
    ]
    meta ""
      pos 279.0
      pos 234.0
  ]
  node [
    id 4
    label "4"
    state "[{}, [aa,b,b]]"
    mutype [
      fixed 1
      hetB 2
    ]
    meta ""
      pos 135.0
      pos 162.0
  ]
  node [
    id 5
    label "5"
    state "[{}, [a,ab,b]]"
    mutype [
      hetA 1
      hetAB 1
      hetB 1
    ]
    meta ""
      pos 63.0
      pos 162.0
  ]
  node [
    id 6
    label "6"
    state "[{}, [a,a,bb]]"
    mutype [
      hetA 2
      fixed 1
    ]
    meta ""
      pos 207.0
      pos 162.0
  ]
  node [
    id 7
    label "7"
    state "[{aa}, [bb]]"
    mutype [
      fixed 2
    ]
    meta ""
      pos 279.0
      pos 162.0
  ]
  node [
    id 8
    label "8"
    state "[{}, [aab,b]]"
    mutype [
      hetB 2
    ]
    meta ""
      pos 99.0
      pos 90.0
  ]
  node [
    id 9
    label "9"
    state "[{}, [aa,bb]]"
    mutype [
      fixed 2
    ]
    meta ""
      pos 243.0
      pos 90.0
  ]
  node [
    id 10
    label "10"
    state "[{}, [ab,ab]]"
    mutype [
      hetAB 2
    ]
    meta ""
      pos 27.0
      pos 90.0
  ]
  node [
    id 11
    label "11"
    state "[{}, [a,abb]]"
    mutype [
      hetA 2
    ]
    meta ""
      pos 171.0
      pos 90.0
  ]
  node [
    id 12
    label "12"
    state "[{}, [aabb]]"
    mutype [
      hetB 2
    ]
    meta "sink"
      pos 135.0
      pos 18.0
  ]
  edge [
    source 0
    target 1
    count 1
    event "E"
  ]
  edge [
    source 0
    target 2
    count 1
    event "C_derived"
  ]
  edge [
    source 0
    target 3
    count 1
    event "C_ancestor"
  ]
  edge [
    source 1
    target 4
    count 1
    event "C_ancestor"
  ]
  edge [
    source 1
    target 5
    count 4
    event "C_ancestor"
  ]
  edge [
    source 1
    target 6
    count 1
    event "C_ancestor"
  ]
  edge [
    source 2
    target 4
    count 1
    event "E"
  ]
  edge [
    source 2
    target 7
    count 1
    event "C_ancestor"
  ]
  edge [
    source 3
    target 6
    count 1
    event "E"
  ]
  edge [
    source 3
    target 7
    count 1
    event "C_derived"
  ]
  edge [
    source 4
    target 8
    count 2
    event "C_ancestor"
  ]
  edge [
    source 4
    target 9
    count 1
    event "C_ancestor"
  ]
  edge [
    source 5
    target 8
    count 1
    event "C_ancestor"
  ]
  edge [
    source 5
    target 10
    count 1
    event "C_ancestor"
  ]
  edge [
    source 5
    target 11
    count 1
    event "C_ancestor"
  ]
  edge [
    source 6
    target 9
    count 1
    event "C_ancestor"
  ]
  edge [
    source 6
    target 11
    count 2
    event "C_ancestor"
  ]
  edge [
    source 7
    target 9
    count 1
    event "E"
  ]
  edge [
    source 8
    target 12
    count 1
    event "C_ancestor"
  ]
  edge [
    source 9
    target 12
    count 1
    event "C_ancestor"
  ]
  edge [
    source 10
    target 12
    count 1
    event "C_ancestor"
  ]
  edge [
    source 11
    target 12
    count 1
    event "C_ancestor"
  ]
]
