;CFD之道
;创建一个标签页框架，parent为父控件，为panel或其他容器，若为panel则为顶级标签页，若为frame则为子标签页
;label为名字，bool1表示是否创建为标签页，#t表示为标签页，#f表示非标签页
;bool2用于指定是否显示边框，#t为显示，#f为不显示
;(cx-create-frame parent label 'tabbed bool1 'border bool2)

(define apply-cb #f)
(define update-cb #f)

;创建一个名为tabbed-frame的框架，其中放了两个选项卡tab1和tab2
(define bzhan666 (cx-create-panel "tab Demo" apply-cb update-cb))
(define tabbed-frame (cx-create-frame bzhan666 #f 'tabbed #t 'border #t))
(define tab1 (cx-create-frame tabbed-frame "tab1" 'border #t))
(define tab2 (cx-create-frame tabbed-frame "tab2" 'border #t))
(define tab3 (cx-create-frame tabbed-frame "tab3" 'border #t))

;在tab1中放入一个名为table1的table,并在table中放入各种功能控件
(define table1 (cx-create-table tab1 "tab1 panel"))
(define ints (cx-create-integer-entry table1 "int:" 'row 0))
(define reals (cx-create-real-entry table1 "real:" 'row 1))
(define txt (cx-create-text-entry table1 "text:" 'row 2))

;创建一个table2，并在其中放入其他的功能控件
(define table2 (cx-create-table tab2 "tab2 panel"))
(define buttonBox (cx-create-button-box table2 "button box" 'radio-mode #f))
(define checkBox1 (cx-create-toggle-button buttonBox "check box 1"))
(define checkBox2 (cx-create-toggle-button buttonBox "check box 2"))

(cx-set-toggle-button checkBox1 #t)

(cx-show-panel bzhan666)
