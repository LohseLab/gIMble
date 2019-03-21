graph [
  directed 1
  state ""
  meta ""
  node [
    id 0
    label "0"
    state "(('a', 'a'), ('b', 'b'))"
    meta "source"
      pos 376.3
      pos 450.0
  ]
  node [
    id 1
    label "1"
    state "(('a', 'a', 'b', 'b'), ())"
    meta ""
      pos 330.3
      pos 234.0
  ]
  node [
    id 2
    label "2"
    state "(('a',), ('a', 'b', 'b'))"
    meta ""
      pos 256.3
      pos 378.0
  ]
  node [
    id 3
    label "3"
    state "(('aa',), ('b', 'b'))"
    meta ""
      pos 504.3
      pos 306.0
  ]
  node [
    id 4
    label "4"
    state "(('a', 'a'), ('bb',))"
    meta ""
      pos 404.3
      pos 378.0
  ]
  node [
    id 5
    label "5"
    state "(('aa', 'b', 'b'), ())"
    meta ""
      pos 474.3
      pos 162.0
  ]
  node [
    id 6
    label "6"
    state "(('a', 'ab', 'b'), ())"
    meta ""
      pos 114.3
      pos 162.0
  ]
  node [
    id 7
    label "7"
    state "(('a', 'a', 'bb'), ())"
    meta ""
      pos 330.3
      pos 162.0
  ]
  node [
    id 8
    label "8"
    state "((), ('a', 'a', 'b', 'b'))"
    meta ""
      pos 256.3
      pos 306.0
  ]
  node [
    id 9
    label "9"
    state "(('a',), ('ab', 'b'))"
    meta ""
      pos 95.299
      pos 306.0
  ]
  node [
    id 10
    label "10"
    state "(('a',), ('a', 'bb'))"
    meta ""
      pos 184.3
      pos 306.0
  ]
  node [
    id 11
    label "11"
    state "((), ('aa', 'b', 'b'))"
    meta ""
      pos 440.3
      pos 234.0
  ]
  node [
    id 12
    label "12"
    state "(('aa',), ('bb',))"
    meta ""
      pos 512.3
      pos 234.0
  ]
  node [
    id 13
    label "13"
    state "(('aab', 'b'), ())"
    meta ""
      pos 294.3
      pos 90.0
  ]
  node [
    id 14
    label "14"
    state "(('aa', 'bb'), ())"
    meta ""
      pos 402.3
      pos 90.0
  ]
  node [
    id 15
    label "15"
    state "(('ab', 'ab'), ())"
    meta ""
      pos 78.299
      pos 90.0
  ]
  node [
    id 16
    label "16"
    state "(('a', 'abb'), ())"
    meta ""
      pos 150.3
      pos 90.0
  ]
  node [
    id 17
    label "17"
    state "((), ('a', 'ab', 'b'))"
    meta ""
      pos 148.3
      pos 234.0
  ]
  node [
    id 18
    label "18"
    state "((), ('a', 'a', 'bb'))"
    meta ""
      pos 258.3
      pos 234.0
  ]
  node [
    id 19
    label "19"
    state "(('a',), ('abb',))"
    meta ""
      pos 76.299
      pos 234.0
  ]
  node [
    id 20
    label "20"
    state "((), ('aab', 'b'))"
    meta ""
      pos 258.3
      pos 162.0
  ]
  node [
    id 21
    label "21"
    state "((), ('aa', 'bb'))"
    meta ""
      pos 402.3
      pos 162.0
  ]
  node [
    id 22
    label "22"
    state "(('aabb',), ())"
    meta "sink"
      pos 222.3
      pos 18.0
  ]
  node [
    id 23
    label "23"
    state "((), ('ab', 'ab'))"
    meta ""
      pos 42.299
      pos 162.0
  ]
  node [
    id 24
    label "24"
    state "((), ('a', 'abb'))"
    meta ""
      pos 186.3
      pos 162.0
  ]
  node [
    id 25
    label "25"
    state "((), ('aabb',))"
    meta ""
      pos 222.3
      pos 90.0
  ]
  edge [
    source 0
    target 1
    count 1
    event "R"
  ]
  edge [
    source 0
    target 2
    count 2
    event "M"
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
    target 5
    count 1
    event "{C}"
  ]
  edge [
    source 1
    target 6
    count 4
    event "{C}"
  ]
  edge [
    source 1
    target 7
    count 1
    event "{C}"
  ]
  edge [
    source 2
    target 1
    count 1
    event "R"
  ]
  edge [
    source 2
    target 8
    count 1
    event "M"
  ]
  edge [
    source 2
    target 9
    count 2
    event "[C]"
  ]
  edge [
    source 2
    target 10
    count 1
    event "[C]"
  ]
  edge [
    source 3
    target 5
    count 1
    event "R"
  ]
  edge [
    source 3
    target 11
    count 1
    event "M"
  ]
  edge [
    source 3
    target 12
    count 1
    event "[C]"
  ]
  edge [
    source 4
    target 7
    count 1
    event "R"
  ]
  edge [
    source 4
    target 10
    count 2
    event "M"
  ]
  edge [
    source 4
    target 12
    count 1
    event "{C}"
  ]
  edge [
    source 5
    target 13
    count 2
    event "{C}"
  ]
  edge [
    source 5
    target 14
    count 1
    event "{C}"
  ]
  edge [
    source 6
    target 13
    count 1
    event "{C}"
  ]
  edge [
    source 6
    target 15
    count 1
    event "{C}"
  ]
  edge [
    source 6
    target 16
    count 1
    event "{C}"
  ]
  edge [
    source 7
    target 14
    count 1
    event "{C}"
  ]
  edge [
    source 7
    target 16
    count 2
    event "{C}"
  ]
  edge [
    source 8
    target 1
    count 1
    event "R"
  ]
  edge [
    source 8
    target 11
    count 1
    event "[C]"
  ]
  edge [
    source 8
    target 17
    count 4
    event "[C]"
  ]
  edge [
    source 8
    target 18
    count 1
    event "[C]"
  ]
  edge [
    source 9
    target 6
    count 1
    event "R"
  ]
  edge [
    source 9
    target 17
    count 1
    event "M"
  ]
  edge [
    source 9
    target 19
    count 1
    event "[C]"
  ]
  edge [
    source 10
    target 7
    count 1
    event "R"
  ]
  edge [
    source 10
    target 18
    count 1
    event "M"
  ]
  edge [
    source 10
    target 19
    count 1
    event "[C]"
  ]
  edge [
    source 11
    target 5
    count 1
    event "R"
  ]
  edge [
    source 11
    target 20
    count 2
    event "[C]"
  ]
  edge [
    source 11
    target 21
    count 1
    event "[C]"
  ]
  edge [
    source 12
    target 14
    count 1
    event "R"
  ]
  edge [
    source 12
    target 21
    count 1
    event "M"
  ]
  edge [
    source 13
    target 22
    count 1
    event "{C}"
  ]
  edge [
    source 14
    target 22
    count 1
    event "{C}"
  ]
  edge [
    source 15
    target 22
    count 1
    event "{C}"
  ]
  edge [
    source 16
    target 22
    count 1
    event "{C}"
  ]
  edge [
    source 17
    target 6
    count 1
    event "R"
  ]
  edge [
    source 17
    target 20
    count 1
    event "[C]"
  ]
  edge [
    source 17
    target 23
    count 1
    event "[C]"
  ]
  edge [
    source 17
    target 24
    count 1
    event "[C]"
  ]
  edge [
    source 18
    target 7
    count 1
    event "R"
  ]
  edge [
    source 18
    target 21
    count 1
    event "[C]"
  ]
  edge [
    source 18
    target 24
    count 2
    event "[C]"
  ]
  edge [
    source 19
    target 16
    count 1
    event "R"
  ]
  edge [
    source 19
    target 24
    count 1
    event "M"
  ]
  edge [
    source 20
    target 13
    count 1
    event "R"
  ]
  edge [
    source 20
    target 25
    count 1
    event "[C]"
  ]
  edge [
    source 21
    target 14
    count 1
    event "R"
  ]
  edge [
    source 21
    target 25
    count 1
    event "[C]"
  ]
  edge [
    source 23
    target 15
    count 1
    event "R"
  ]
  edge [
    source 23
    target 25
    count 1
    event "[C]"
  ]
  edge [
    source 24
    target 16
    count 1
    event "R"
  ]
  edge [
    source 24
    target 25
    count 1
    event "[C]"
  ]
  edge [
    source 25
    target 22
    count 1
    event "R"
  ]
]
