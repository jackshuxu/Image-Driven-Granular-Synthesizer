
(princ "LOADING TEST RUNTIME DEBUG VERSION\n")

(load "nyinit")

(defun gran-t (filename &optional
              (duration 1.0)
              (file-pos 0.0)     ;[0:1]        float
              (section-size 0.5);[0:1]         float
              (env-index 1)      ;[1:3]        int
              (env-width 0.125)  ;[0:1]        float
              (density 25)       ;[2:200]      int
              (gran-size 25)     ;[5:200]      int
              (speed 1.0)        ;[-2:2]       float
              (gain 0.0)           ;[-48:6]    float
              (random 0.0)       ;[0:100]      float
              (random-spread 0.0);[0:1]        float
              (pitch 0.0)        ;[-2000:2000] float
              )
      (let ((wave (s-read filename)))
                  (snd-gran wave
                            (snd-read-dur *rslt*)
                            duration
                            file-pos
                            section-size
                            env-index
                            env-width
                            density
                            gran-size
                            speed
                            gain
                            random
                            random-spread
                            pitch
                  )
      )
)

(defun gran-sound (filename duration image-data)
      (mult (gran-t filename duration 0.25 0.28 2 0.1 200 15 1.8 0.0 88 0.5 202.0) (pwl 0.05 0.9 2.45 0.7 duration 0.0))
)


;; (snd-juce-import)
;; ;; (play (gran_t "Piano_C4.wav" 10.0 0.1 0.5 1 0.125 50 25 -0.5 0.0 50.0 0.5 500.0))
;; ;; (timed-seq (seq (gran-t "Piano_C4.wav" 2.0 0.1 0.5 1 0.125 50 25 -0.5 0.0 50.0 0.5 500.0)))
;; ;;            2.0 4.0 (gran-t "Piano_C4.wav" 2.0 0.1 0.5 1 0.125 50 25 -0.5 0.0 50.0 0.5 -500.0))
      
;; ;; )

;; (play (timed-seq '((0.0 2.0 (gran-t "Piano_C4.wav" 2.0 0.1 0.5 1 0.125 50 25 -0.5 0.0 50.0 0.5 500.0))
;;                          (2.0 4.0 (gran-t "Piano_C4.wav" 2.0 0.1 0.5 1 0.125 50 25 -0.5 0.0 50.0 0.5 -500.0))
;;                         )
;;             )
;; )

;; (snd-test (osc c4))