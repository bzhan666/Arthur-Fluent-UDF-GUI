(define apply-cb #t)
(define update-cb #f)

(define bzhan666 (cx-create-panel "bzhan666" apply-cb update-cb))
(define table)
(define checkbox1)
(define checkbox2)
(define checkbox3)
(define checkbox4)
;(define isBool)

(set! table(cx-create-table bzhan666 "cfd"))

(define buttonBox(cx-create-button-box table "Button box" 'radio-mode #t))
(set! checkbox1 (cx-create-toggle-button buttonBox "Pengyuyan1"))
(set! checkbox2 (cx-create-toggle-button buttonBox "Pengyuyan2"))
(set! checkbox3 (cx-create-toggle-button buttonBox "Pengyuyan3"))
(set! checkbox4 (cx-create-toggle-button buttonBox "Pengyuyan4"))

;(cx-set-toggle-button checkbox3 #t)

(cx-show-panel bzhan666)