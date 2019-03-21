graph [
  directed 1
  state ""
  meta ""
  node [
    id 0
    label "0"
    state "(('a', 'a'), ('b', 'b'))"
    meta "source"
      pos 198.0
      pos 378.0
  ]
  node [
    id 1
    label "1"
    state "(('a',), ('a', 'b', 'b'))"
    meta ""
      pos 107.0
      pos 306.0
  ]
  node [
    id 2
    label "2"
    state "((), ('a', 'a', 'b', 'b'))"
    meta ""
      pos 171.0
      pos 234.0
  ]
  node [
    id 3
    label "3"
    state "(('aa',), ('b', 'b'))"
    meta ""
      pos 298.0
      pos 306.0
  ]
  node [
    id 4
    label "4"
    state "(('a', 'a'), ('bb',))"
    meta ""
      pos 226.0
      pos 306.0
  ]
  node [
    id 5
    label "5"
    state "(('a',), ('ab', 'b'))"
    meta ""
      pos 27.0
      pos 234.0
  ]
  node [
    id 6
    label "6"
    state "(('a',), ('a', 'bb'))"
    meta ""
      pos 99.0
      pos 234.0
  ]
  node [
    id 7
    label "7"
    state "((), ('aa', 'b', 'b'))"
    meta ""
      pos 253.0
      pos 162.0
  ]
  node [
    id 8
    label "8"
    state "((), ('a', 'ab', 'b'))"
    meta ""
      pos 109.0
      pos 162.0
  ]
  node [
    id 9
    label "9"
    state "((), ('a', 'a', 'bb'))"
    meta ""
      pos 181.0
      pos 162.0
  ]
  node [
    id 10
    label "10"
    state "(('aa',), ('bb',))"
    meta ""
      pos 319.0
      pos 234.0
  ]
  node [
    id 11
    label "11"
    state "(('a',), ('abb',))"
    meta ""
      pos 37.0
      pos 162.0
  ]
  node [
    id 12
    label "12"
    state "((), ('aab', 'b'))"
    meta ""
      pos 181.0
      pos 90.0
  ]
  node [
    id 13
    label "13"
    state "((), ('aa', 'bb'))"
    meta ""
      pos 253.0
      pos 90.0
  ]
  node [
    id 14
    label "14"
    state "((), ('ab', 'ab'))"
    meta ""
      pos 37.0
      pos 90.0
  ]
  node [
    id 15
    label "15"
    state "((), ('a', 'abb'))"
    meta ""
      pos 109.0
      pos 90.0
  ]
  node [
    id 16
    label "16"
    state "((), ('aabb',))"
    meta "sink"
      pos 145.0
      pos 18.0
  ]
  edge [
    source 0
    target 1
    count 2
    event "M"
  ]
  edge [
    source 0
    target 2
    count 1
    event "E"
  ]
  edge [
    source 0
    target 3
    count 1
    event "{C}"
  ]
  edge [
    source 0
    target 4
    count 1
    event "[C]"
  ]
  edge [
    source 1
    target 2
    count 1
    event "E"
  ]
  edge [
    source 1
    target 5
    count 2
    event "[C]"
  ]
  edge [
    source 1
    target 6
    count 1
    event "[C]"
  ]
  edge [
    source 2
    target 7
    count 1
    event "[C]"
  ]
  edge [
    source 2
    target 8
    count 4
    event "[C]"
  ]
  edge [
    source 2
    target 9
    count 1
    event "[C]"
  ]
  edge [
    source 3
    target 7
    count 1
    event "E"
  ]
  edge [
    source 3
    target 10
    count 1
    event "[C]"
  ]
  edge [
    source 4
    target 6
    count 2
    event "M"
  ]
  edge [
    source 4
    target 9
    count 1
    event "E"
  ]
  edge [
    source 4
    target 10
    count 1
    event "{C}"
  ]
  edge [
    source 5
    target 8
    count 1
    event "E"
  ]
  edge [
    source 5
    target 11
    count 1
    event "[C]"
  ]
  edge [
    source 6
    target 9
    count 1
    event "E"
  ]
  edge [
    source 6
    target 11
    count 1
    event "[C]"
  ]
  edge [
    source 7
    target 12
    count 2
    event "[C]"
  ]
  edge [
    source 7
    target 13
    count 1
    event "[C]"
  ]
  edge [
    source 8
    target 12
    count 1
    event "[C]"
  ]
  edge [
    source 8
    target 14
    count 1
    event "[C]"
  ]
  edge [
    source 8
    target 15
    count 1
    event "[C]"
  ]
  edge [
    source 9
    target 13
    count 1
    event "[C]"
  ]
  edge [
    source 9
    target 15
    count 2
    event "[C]"
  ]
  edge [
    source 10
    target 13
    count 1
    event "E"
  ]
  edge [
    source 11
    target 15
    count 1
    event "E"
  ]
  edge [
    source 12
    target 16
    count 1
    event "[C]"
  ]
  edge [
    source 13
    target 16
    count 1
    event "[C]"
  ]
  edge [
    source 14
    target 16
    count 1
    event "[C]"
  ]
  edge [
    source 15
    target 16
    count 1
    event "[C]"
  ]
]
