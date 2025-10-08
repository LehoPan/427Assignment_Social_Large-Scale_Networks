graph [
  directed 0
  node [
    id 0
    label "A"
    color "red"
  ]
  node [
    id 1
    label "B"
    color "red"
  ]
  node [
    id 2
    label "C"
    color "blue"
  ]
  node [
    id 3
    label "D"
    color "blue"
  ]
  node [
    id 4
    label "E"
    color "red"
  ]
  edge [
    source 0
    target 1
    sign 1
  ]
  edge [
    source 1
    target 2
    sign -1
  ]
  edge [
    source 2
    target 3
    sign 1
  ]
  edge [
    source 3
    target 0
    sign -1
  ]
  edge [
    source 0
    target 4
    sign 1
  ]
  edge [
    source 4
    target 1
    sign 1
  ]
]
