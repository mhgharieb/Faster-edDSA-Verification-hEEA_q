BINGCD=bingcd/src
EEA_q=EEA_q
# SAFEGCD=inverse25519skylake-20201231
SAFEGCD=inverse25519skylake-20210110

.PHONY: all bingcd EEA_q safegcd

all: bingcd EEA_q safegcd

bingcd:
	$(MAKE) -C $(BINGCD)

EEA_q:
	$(MAKE) -C $(EEA_q)

safegcd:
	$(MAKE) -C $(SAFEGCD)

clean:
	$(MAKE) -C $(BINGCD) clean && \
	$(MAKE) -C $(EEA_q)  clean && \
	$(MAKE) -C $(SAFEGCD) clean