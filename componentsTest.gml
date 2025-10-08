graph [
  directed 0

  # Cluster 1
  node [
    id 0
    label "A"
  ]
  node [
    id 1
    label "B"
  ]
  node [
    id 2
    label "C"
  ]

  edge [
    source 0
    target 1
  ]
  edge [
    source 1
    target 2
  ]
  edge [
    source 0
    target 2
  ]

  # Cluster 2
  node [
    id 3
    label "D"
  ]
  node [
    id 4
    label "E"
  ]
  node [
    id 5
    label "F"
  ]

  edge [
    source 3
    target 4
  ]
  edge [
    source 4
    target 5
  ]
  edge [
    source 3
    target 5
  ]

  # Sparse inter-cluster connection
  edge [
    source 2
    target 3
  ]
]
