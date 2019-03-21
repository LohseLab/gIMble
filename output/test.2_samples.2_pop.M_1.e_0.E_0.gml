graph [
  directed 1
  state ""
  meta ""
  node [
    id 0
    label "0"
    state "(('a', 'a'), ('b', 'b'))"
    meta "source"
      pos 111.48
      pos 378.0
  ]
  node [
    id 1
    label "1"
    state "(('a',), ('a', 'b', 'b'))"
    meta ""
      pos 183.48
      pos 306.0
  ]
  node [
    id 2
    label "2"
    state "(('aa',), ('b', 'b'))"
    meta ""
      pos 39.477
      pos 306.0
  ]
  node [
    id 3
    label "3"
    state "(('a', 'a'), ('bb',))"
    meta ""
      pos 111.48
      pos 306.0
  ]
  node [
    id 4
    label "4"
    state "((), ('a', 'a', 'b', 'b'))"
    meta ""
      pos 111.48
      pos 234.0
  ]
  node [
    id 5
    label "5"
    state "(('a',), ('ab', 'b'))"
    meta ""
      pos 255.48
      pos 234.0
  ]
  node [
    id 6
    label "6"
    state "(('a',), ('a', 'bb'))"
    meta ""
      pos 183.48
      pos 234.0
  ]
  node [
    id 7
    label "7"
    state "((), ('aa', 'b', 'b'))"
    meta ""
      pos 39.477
      pos 162.0
  ]
  node [
    id 8
    label "8"
    state "(('aa',), ('bb',))"
    meta ""
      pos 39.477
      pos 234.0
  ]
  node [
    id 9
    label "9"
    state "((), ('a', 'ab', 'b'))"
    meta ""
      pos 183.48
      pos 162.0
  ]
  node [
    id 10
    label "10"
    state "((), ('a', 'a', 'bb'))"
    meta ""
      pos 111.48
      pos 162.0
  ]
  node [
    id 11
    label "11"
    state "(('a',), ('abb',))"
    meta ""
      pos 255.48
      pos 162.0
  ]
  node [
    id 12
    label "12"
    state "((), ('aab', 'b'))"
    meta ""
      pos 111.48
      pos 90.0
  ]
  node [
    id 13
    label "13"
    state "((), ('aa', 'bb'))"
    meta ""
      pos 39.477
      pos 90.0
  ]
  node [
    id 14
    label "14"
    state "((), ('ab', 'ab'))"
    meta ""
      pos 183.48
      pos 90.0
  ]
  node [
    id 15
    label "15"
    state "((), ('a', 'abb'))"
    meta ""
      pos 255.48
      pos 90.0
  ]
  node [
    id 16
    label "16"
    state "((), ('aabb',))"
    meta "sink"
      pos 147.48
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
    event "{C}"
  ]
  edge [
    source 0
    target 3
    count 1
    event "[C]"
  ]
  edge [
    source 1
    target 4
    count 1
    event "M"
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
    event "M"
  ]
  edge [
    source 2
    target 8
    count 1
    event "[C]"
  ]
  edge [
    source 3
    target 6
    count 2
    event "M"
  ]
  edge [
    source 3
    target 8
    count 1
    event "{C}"
  ]
  edge [
    source 4
    target 7
    count 1
    event "[C]"
  ]
  edge [
    source 4
    target 9
    count 4
    event "[C]"
  ]
  edge [
    source 4
    target 10
    count 1
    event "[C]"
  ]
  edge [
    source 5
    target 9
    count 1
    event "M"
  ]
  edge [
    source 5
    target 11
    count 1
    event "[C]"
  ]
  edge [
    source 6
    target 10
    count 1
    event "M"
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
    target 13
    count 1
    event "M"
  ]
  edge [
    source 9
    target 12
    count 1
    event "[C]"
  ]
  edge [
    source 9
    target 14
    count 1
    event "[C]"
  ]
  edge [
    source 9
    target 15
    count 1
    event "[C]"
  ]
  edge [
    source 10
    target 13
    count 1
    event "[C]"
  ]
  edge [
    source 10
    target 15
    count 2
    event "[C]"
  ]
  edge [
    source 11
    target 15
    count 1
    event "M"
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
