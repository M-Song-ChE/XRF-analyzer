# XRF Data Analyzer — 개발 현황

GitHub: https://github.com/M-Song-ChE/XRF-analyzer  
최종 업데이트: 2026-05-17

---

## 개요

XRF(X-ray Fluorescence) 스팟맵 CSV 데이터를 불러와 mass% → atomic fraction으로 변환하고,  
선택한 원소들끼리의 평균 조성비와 표준편차를 계산해 보여주는 tkinter GUI 앱.

---

## 구현 완료 기능

### 데이터 처리
- CSV 파일 다중 로드 (파일마다 독립적으로 관리)
- CSV 포맷 자동 인식 — X, Y 헤더 행 위치를 직접 탐색하므로 앞부분 메타데이터 행 수에 무관
- mass% → atomic fraction 변환 (스팟별 moles = mass% / 원자량, 선택 원소끼리 정규화)
- 섹션별 조성비 계산:
  - 선택한 원소들의 moles 합계(subtotal)가 **정확히 0.0**인 섹션만 "Undefined"로 제외
  - 하나라도 > 0이면 포함 (예: Pt=0 & Ni>0 → 포함)
- 평균(mean) 및 표본 표준편차(std, ddof=1) 계산

### UI — 주기율표 패널
- 전체 118개 원소 인터랙티브 주기율표 표시
- 원소 색상 — 카테고리별 구분 (alkali, alkaline-earth, transition, post-transition, metalloid, nonmetal, halogen, noble-gas, lanthanide, actinide)
- 3가지 버튼 상태: 포함(vivid), 검출됨-미포함(dim+groove), 미검출(flat/비활성)
- 클릭으로 토글 — 선택 원소끼리만 재정규화해 조성 계산
- **줌 인/아웃**: 툴바 `−` / `+` 버튼 및 마우스 휠 (범위: 자연 크기의 30%–300%)
- 줌 시 주기율표가 수평 중앙 기준으로 확대/축소
- 창 크기 변경 시 자동 리사이즈 (zoom factor 적용)
- 기본 줌: 자연 크기의 50%

### UI — 원소별 통계 표 (중간 패널)
- 검출된 원소 전체 표시
- 컬럼: Element / Z / Mean at% (renorm) / ±σ / N spots / Mean mass%
- 포함 원소: 녹색 배경 + 볼드, 미포함 원소: 회색 + 괄호로 참고값 표시
- 컬럼 헤더 클릭으로 정렬

### UI — 퍼파일 조성 표 (하단 패널)
- 파일마다 한 행 — 파일명, 조성(합금 표기식), Total/Defined/Undefined 섹션 수, 원소별 at% ± σ
- 합금 표기식: `Pt₄₅.₃ Ni₅₄.₇` 형식 (Unicode 아래 첨자)
- **CSV 내보내기** 버튼
- 파일 클릭 시 해당 파일 단독 뷰 / "View All Files Combined"로 전체 보기

### UI — 레이아웃
- 좌측 사이드바 (파일 목록, 원소 선택 컨트롤, 범례) / 우측 3-패널 수직 분할
- `ttk.PanedWindow`로 각 패널 크기 드래그 조절 가능
- 앱 시작 시 사시 위치 자동 설정 (주기율표 ~1/4, 원소표 ~55%, 퍼파일표 나머지)

---

## 파일 구조

```
XRF data extractor/
├── xrf_analyzer.py      # 메인 단일 파일 앱 (~1000줄)
├── DEV_STATUS.md        # 이 파일
├── PtNi_Rec_2.csv       # 샘플 데이터
├── .gitignore
└── .venv/               # Python 가상환경 (numpy, tkinter)
```

---

## 계산 방법 요약

```
스팟 하나당:
  moles_i = mass%_i / atomic_mass_i   (선택 원소만)
  subtotal = sum(moles_i)

  if subtotal == 0.0:  → Undefined (제외)
  else:  at%_i = moles_i / subtotal

파일 전체:
  mean_at%_i = mean(at%_i across valid spots)
  σ_i        = std(at%_i across valid spots, ddof=1)
```

---

## 커밋 이력

| 해시 | 내용 |
|------|------|
| `dd5a438` | 합금 표기식 박스 제거, 원소별 표 복원, PT 기본 크기 절반 |
| `99fd439` | 레이아웃 수정: PT 높이 고정, PT 중앙 정렬, 퍼파일 컬럼 stretch |
| `566f81e` | PT 줌 컨트롤 추가, 테이블 폰트 조정, 컬럼 너비 균등화 |
| `fddcd55` | 초기 커밋 |

---

## 알려진 이슈 / 향후 개선 가능 사항

- [ ] 퍼파일 조성 표에서 조성 열 텍스트가 길어질 경우 잘릴 수 있음 (tooltip 고려)
- [ ] 파일 제거 후 포커스 처리 개선 여지
- [ ] 원소 통계 표 정렬 컬럼명이 내부 컬럼 ID와 불일치 (minor)
- [ ] 향후: 스팟 좌표 기반 히트맵 시각화
- [ ] 향후: 여러 파일 간 조성 비교 그래프
