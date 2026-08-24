* Integrality declared the two ways MPS allows: a MARKER block around a run of
* columns, and a BOUNDS entry naming a column directly. Neither survives the
* parser, so the reader recovers both from the source text.
NAME          markers
ROWS
 N  COST
 G  C1
 L  C2
COLUMNS
    MARKER                 'MARKER'                 'INTORG'
    XINT      COST         3.0        C1           1.0
    XINT      C2           1.0
    MARKER                 'MARKER'                 'INTEND'
    XCONT     COST         1.0        C1           1.0
    XCONT     C2           1.0
    XBIN      COST         7.0        C1           2.0
RHS
    RHS       C1           2.5        C2           9.0
BOUNDS
 UI BND       XINT         9.0
 UP BND       XCONT        4.0
 BV BND       XBIN
ENDATA
